/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "polyLangevin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyLangevin, 0);

addToRunTimeSelectionTable(polyStateController, polyLangevin, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyLangevin::polyLangevin
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    molIds_(),
    storedVelocities_(),
    gamma_(readScalar(propsDict_.lookup("gamma"))),
    controlOn_(true),    
    switchOffAfterTime_(false),
    currentTime_(0.0),
    switchOffTime_(GREAT)    
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

//     singleValueController() = true;
    
    deltaT_ = time_.deltaT().value(); 
    
    // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    if (propsDict_.found("temperature"))
    {
        temperature_ = readScalar(propsDict_.lookup("temperature"));
    }
    else
    {
        const reducedUnits& rU = molCloud_.redUnits();
        
        temperature_ = readScalar(propsDict_.lookup("temperatureKelvin"));
        
        temperature_ /= rU.refTemp();
    }
//     friction_ = 1.0 - (0.5*gamma_*deltaT_);
//     noise_ = sqrt((6.0*gamma_*temperature_) / deltaT_);
    
//     Info<< "friction coefficient: " << friction_
//         << ", noise: " << noise_
//         << endl;    
    sigma1_ = (1.0 - exp(-1.0*gamma_*deltaT_))/(1.0*gamma_);
    sigma2_ = (1.0 - exp(-2.0*gamma_*deltaT_))/(2.0*gamma_);    
    
    if (propsDict_.found("switchOffAfterTime"))
    {   
        switchOffAfterTime_ = true;
        switchOffTime_ = readScalar(propsDict_.lookup("switchOffAfterTime"));
    }    
    
    elapsedtime_ = 0.0;
    
    variableGamma_ = false;
    
    if (propsDict_.found("variableGamma"))
    {   
        variableGamma_ = Switch(propsDict_.lookup("variableGamma"));
        
        if(variableGamma_)
        {
            deltaGamma_ = readScalar(propsDict_.lookup("deltaGamma"));            
            timeInterval_ = readScalar(propsDict_.lookup("timeInterval"));            
        }
    }    
//     storedVelocities_.setSize(molCloud_.size(), vector::zero);
//     trackingNumbers_.setSize(molCloud_.size(), -1);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyLangevin::~polyLangevin()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyLangevin::initialConfiguration()
{
    label cloudSize = molCloud_.moleculeTracking().getMaxTrackingNumber();
    trackingNumbers_.setSize(cloudSize);
    storedVelocities_.setSize(cloudSize, vector::zero);
}

void polyLangevin::controlBeforeVelocityI()
{}

void polyLangevin::measureTemperature()
{
 
    scalar kE = 0.0;
    scalar dof = 0.0;
    scalar angularKeSum = 0.0;

    // - kinetic energy (from thermal velocities)
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)            
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());

                kE += 0.5*massI*magSqr(mol().v());
                dof += molCloud_.cP().degreesOfFreedom(mol().id());

                const diagTensor& molMoI(molCloud_.cP().momentOfInertia(mol().id()));

                // angular speed 
                const vector& molOmega(inv(molMoI) & mol().pi());
                angularKeSum += 0.5*(molOmega & molMoI & molOmega);
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(kE, sumOp<scalar>());
        reduce(dof, sumOp<scalar>());
        reduce(angularKeSum, sumOp<scalar>());
    }

    scalar tempMeasI = 0.0;

    if(dof > 0.0)
    {
        const scalar& kB = molCloud_.redUnits().kB();

        tempMeasI = (2.0*(kE+angularKeSum))/(kB*dof);

        const reducedUnits& rU = molCloud_.redUnits();

        Info << " (Langevin): T = " << tempMeasI << " (reduced units) "
            << " T = " << tempMeasI*rU.refTemp() << " (SI units) "
            << endl;
    }

}

void polyLangevin::controlBeforeMove()
{
    if(controlOn_)
    {
        storedVelocities_ = vector::zero;
        
//         Info << "temperatureLangevin: control 1" << endl;
        
        scalar sigma1 = sigma1_;
        scalar sigma2 = sigma2_;
        
        const scalar& kB = molCloud_.redUnits().kB();
//         label i = 0;

        DynamicList<label> TNs;
        DynamicList<vector> velocities;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            
            if(!mol().frozen())
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());
                
                scalar termA = sqrt(2*gamma_*massI*kB*temperature_)/massI;
                
                vector z1 = vector::zero;
                vector z2 = vector::zero;   
                
                // assumption
    //             scalar sigma1 = (1.0 - exp(-1.0*gamma_*deltaT_))/(1.0*gamma_);
    //             scalar sigma2 = (1.0 - exp(-2.0*gamma_*deltaT_))/(2.0*gamma_);
                
                vector R1 = vector
                            (
                                molCloud_.rndGen().GaussNormalMD<scalar>(),
                                molCloud_.rndGen().GaussNormalMD<scalar>(),
                                molCloud_.rndGen().GaussNormalMD<scalar>()
                            );
                            
                vector R2 = vector
                            (
                                molCloud_.rndGen().GaussNormalMD<scalar>(),
                                molCloud_.rndGen().GaussNormalMD<scalar>(),
                                molCloud_.rndGen().GaussNormalMD<scalar>()
                            );
                            
                z1 = sqrt(sigma2)*R1;
                z2 = R1*(sigma1-sigma2)/sqrt(sigma2) + R2*sqrt(deltaT_ - (sigma1*sigma1/sigma2));
                
                
                const label& tN = mol().trackingNumber();
                
                storedVelocities_[tN] = exp(-gamma_*deltaT_)*mol().v() + termA*z1;
                
                TNs.append(tN);
                velocities.append(storedVelocities_[tN]);
                
                // New velocity to temporary hack the move function
                vector newVel = (
                                    ( ( 1.0 - exp(-gamma_*deltaT_) )*mol().v()/gamma_ )
                                    + (termA*z2/gamma_)
                                    
                                )/deltaT_;
                                    
                mol().v() = newVel;
            }
        }
        
        TNs.shrink();
        velocities.shrink();
        
        if(Pstream::parRun())
        {        
            List<label> tNsTransf(TNs.size());        
            List<vector> velTransf(TNs.size());        
            tNsTransf.transfer(TNs);
            velTransf.transfer(velocities);
            
            //- sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << tNsTransf << velTransf;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    List<label> tNsProc;
                    List<vector> velProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> tNsProc >> velProc;
                    }
            
                    forAll(tNsProc, i)
                    {
                        storedVelocities_[tNsProc[i]] = velProc[i];
                    }
                }
            }
        }
//         if(Pstream::parRun())
//         {
//             forAll(storedVelocities_, i)
//             {
//                 reduce(storedVelocities_[i], sumOp<vector>());
//             }
//         }
    }
    else
    {
//         Info << "temperatureLangevin: OFF" << endl;
    }
}

void polyLangevin::controlBeforeForces()
{}

void polyLangevin::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyLangevin::controlAfterForces()
{}


void polyLangevin::controlAfterVelocityII()
{
    measureTemperature();
    
    if(controlOn_)
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            const label& tN = mol().trackingNumber();
            
            if(!mol().frozen())
            {
                mol().v() = storedVelocities_[tN] + 0.5*deltaT_*mol().a();
            }
        }
    }
    
    
    if(switchOffAfterTime_)
    {
        currentTime_ += deltaT_;

        if(currentTime_ >= switchOffTime_)
        {
            controlOn_ = false;
        }
    }
    
    
    if(variableGamma_)
    {
        elapsedtime_ += deltaT_;
        
        if(elapsedtime_ >= timeInterval_)
        {
            elapsedtime_ = 0.0;
            gamma_ += deltaGamma_;
        }        
        
//         Info << "gamma = " << gamma_ << endl;
    }
    
}



void polyLangevin::calculateProperties()
{}

void polyLangevin::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}



void polyLangevin::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
}




} // End namespace Foam

// ************************************************************************* //
