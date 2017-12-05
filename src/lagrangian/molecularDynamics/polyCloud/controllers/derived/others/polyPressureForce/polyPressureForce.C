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

#include "polyPressureForce.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyPressureForce, 0);

addToRunTimeSelectionTable(polyStateController, polyPressureForce, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPressureForce::polyPressureForce
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fileName_(propsDict_.lookup("fileName")),    
    molIds_()
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    
    d_ = propsDict_.lookup("forceDirection");
    d_ /= mag(d_);
    
    area_ = readScalar(propsDict_.lookup("area"));

    targetPressure_ = readScalar(propsDict_.lookup("pressureMPa"));    
    
    targetPressure_ *= 1e6;
    
    targetPressure_ /= molCloud_.redUnits().refPressure();
    
    Info << "target Pressure RU " <<  targetPressure_ << endl;
    
    readProperties();
    
    endTime_ = 0;
    
    if(propsDict_.found("endTime"))
    {
        endTime_ = readScalar(propsDict_.lookup("endTime"));
    }
    
    Info<< "polyPressureForce: max force = " << targetPressure_*area_*d_ << endl;
 
    if(endTime_ > 0)
    {
        force_ = vector::zero;
        dt_ = (targetPressure_*area_)/endTime_;
    }
    else 
    {
        endTime_ = 0.0;
        force_ = targetPressure_*area_*d_;  
    }

    currentTime_ = 0.0;
    
    deltaT_ = t.deltaT().value();
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPressureForce::~polyPressureForce()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyPressureForce::initialConfiguration()
{
    DynamicList<label> trackingNumbers;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                trackingNumbers.append(mol().trackingNumber());
                
                // change from non-frozen to frozen 
                
                mol().special() = 0;
            }
        }
    } 
    
    trackingNumbers.shrink();
    
    // parallel communication

    if(Pstream::parRun())
    { 
         List<label> tNs (trackingNumbers.size(), 0);
        
        forAll(trackingNumbers, i)
        {
            tNs[i] = trackingNumbers[i];
        }
        
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << tNs;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<label> tNsProc;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> tNsProc;
                }
                
                forAll(tNsProc, i)
                {
                    trackingNumbers.append(tNsProc[i]);
                }
            }
        }    
    }
    
    trackingNumbers.shrink();
    
    trackingNumbers_.setSize(trackingNumbers.size());
    
    trackingNumbers_.transfer(trackingNumbers);
    
    centreOfMass();
    
    // set initial velocities on all carbon atoms the same
    setVelocities();
    
    controlAfterForces();
}

void polyPressureForce::setVelocities()
{
    vector velocity(vector::zero);
    scalar nMols = 0.0;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                velocity += mol().v();
                nMols += 1.0;
            }
        }
    }
    
    // parallel communication

    if(Pstream::parRun())
    { 
        reduce(nMols, sumOp<scalar>());
        reduce(velocity, sumOp<vector>());        
    }
    
    if(nMols > 0)
    {
        velocity /= nMols;
    }
    else
    {
        Info << "WARNING: polyPressureForce: no mols - something is wrong " << endl;
    }
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                mol().v() = (velocity & d_)*d_;
            }
        }   
    }
    
    velocity_ = velocity;
}

void polyPressureForce::centreOfMass()
{
    vector centreOfMass(vector::zero);
    scalar mass = 0.0;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());                
                
                centreOfMass += mol().position()*massI;
                mass += massI;
            }
        }
    } 
 
    // parallel communication

    if(Pstream::parRun())
    { 
        reduce(mass, sumOp<scalar>());
        reduce(centreOfMass, sumOp<vector>());        
    }
    
    position_ = centreOfMass/mass;
    
    mass_ = mass;
}


void polyPressureForce::controlBeforeVelocityI()
{}

void polyPressureForce::controlBeforeMove()
{}

void polyPressureForce::controlBeforeForces()
{}

void polyPressureForce::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyPressureForce::controlAfterForces()
{
    centreOfMass();
    
    setVelocities();
    
    currentTime_ += deltaT_;

    if(currentTime_ < endTime_)
    {
        force_ = dt_*d_*currentTime_;
    }
    else
    {
        force_ = targetPressure_*area_*d_;
    }
    
    Info<< "polyPressureForce: force = " << force_ 
        << ", velocity = " << velocity_
        << ", acc = " << netAcc_
        << ", position = " << position_
        << endl;

    // Get the total force on the object after intermolecular force computation
    
    {
        netForce_ = vector::zero;
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                forAll(mol().siteForces(), i)
                {
                    netForce_ += mol().siteForces()[i];
                }
            }
        }
        
        // parallel communication        
        if(Pstream::parRun())
        { 
            reduce(netForce_, sumOp<vector>());        
        }        
        
        // plus external pressure force
        netForce_ += force_;
        
        netAcc_ = netForce_/mass_;
    }       
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                forAll(mol().siteForces(), i)
                {
                    mol().siteForces() = netForce_;
                }

                mol().a() = (netAcc_ & d_)*d_;
                mol().pi() = vector::zero;
                mol().tau() = vector::zero;
            }
        }      
    }    
}

void polyPressureForce::controlAfterVelocityII()
{}


void polyPressureForce::calculateProperties()
{}


void polyPressureForce::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, time_.time().timeOutputValue());
            vectorField positions(1, position_);
            vectorField netForce(1, netForce_);            
            vectorField velocity(1, velocity_);     
            
            writeTimeData
            (
                fixedPathName,
                "polyPressureForce_"+fileName_+"_position.xyz",
                timeField,
                positions,
                true
            );            
            
            writeTimeData
            (
                fixedPathName,
                "polyPressureForce_"+fileName_+"_netForce.xyz",
                timeField,
                netForce,
                true
            );
            
            writeTimeData
            (
                fixedPathName,
                "polyPressureForce_"+fileName_+"_velocity.xyz",
                timeField,
                velocity,
                true
            );            
        }
    }
}

void polyPressureForce::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    readProperties();
}


void polyPressureForce::readProperties()
{}




} // End namespace Foam

// ************************************************************************* //
