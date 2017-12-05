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

#include "polyTemperatureBerendsenBinsNew.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyTemperatureBerendsenBinsNew, 0);

addToRunTimeSelectionTable(polyStateController, polyTemperatureBerendsenBinsNew, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyTemperatureBerendsenBinsNew::polyTemperatureBerendsenBinsNew
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, /*mesh,*/ molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
 	tauT_(readScalar(propsDict_.lookup("tauT"))),
    componentControl_(false),
    X_(false),
    Y_(false),
    Z_(false),
    peculiar_(false),

    molIds_(),

    averagingCounter_(0.0),
    resetFieldsAtOutput_(true),
    output_(false)
{
    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;


    temperature_ = readScalar(propsDict_.lookup("temperature"));

    if (propsDict_.found("componentControl"))
    {
        componentControl_ = Switch(propsDict_.lookup("componentControl"));

        if(componentControl_)
        {
            if (propsDict_.found("X"))
            {
                X_ = Switch(propsDict_.lookup("X"));
            }
    
            if (propsDict_.found("Y"))
            {
                Y_ = Switch(propsDict_.lookup("Y"));
            }
    
            if (propsDict_.found("Z"))
            {
                Z_ = Switch(propsDict_.lookup("Z"));
            }

            if(!X_ && !Y_ && !Z_)
            {
                FatalErrorIn("polyTemperatureBerendsenBinsNew::polyTemperatureBerendsenBinsNew()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }
    
    measureFullTemperature_ = false;
    
    if (propsDict_.found("measureFullTemperature"))
    {
        measureFullTemperature_ = Switch(propsDict_.lookup("measureFullTemperature"));
    }    

    if (propsDict_.found("peculiar"))
    {
        peculiar_ = Switch(propsDict_.lookup("peculiar"));
    }

    if (propsDict_.found("resetFieldsAtOutput"))
    {
        resetFieldsAtOutput_ = Switch(propsDict_.lookup("resetFieldsAtOutput"));
    }
    
    angularControl_ = false;
    
    if (propsDict_.found("angularControl"))
    {
        angularControl_ = Switch(propsDict_.lookup("angularControl"));
    }

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    binModel_ =  autoPtr<binModel>
    (
        binModel::New(mesh_, propsDict_)
    );


    const label& nBins = binModel_->nBins();
    mass_.setSize(nBins, 0.0);
    momentum_.setSize(nBins, vector::zero);
    velocity_.setSize(nBins, vector::zero);
    kE_.setSize(nBins, 0.0);
    angularKe_.setSize(nBins, 0.0);
    dof_.setSize(nBins, 0.0);
    measuredT_.setSize(nBins, temperature_);
    chi_.setSize(nBins, 1.0);
    chiLin_.setSize(nBins, 1.0);
    chiAng_.setSize(nBins, 1.0);
    
    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTemperatureBerendsenBinsNew::~polyTemperatureBerendsenBinsNew()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyTemperatureBerendsenBinsNew::initialConfiguration()
{
    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];

    scalarField mass(mass_.size(), 0.0);
    vectorField momentum(mass_.size(), vector::zero);

    forAll(cells, c)
    {
        const label& cell = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cell];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cell);

            if(n != -1)
            {
                if((findIndex(molIds_, molI->id()) != -1) || measureFullTemperature_)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

//                     Info << "molI->position() "  << molI->position() 
//                     << ", n = " << n
//                     << endl;
                    
                    mass[n] += massI;
                    momentum[n] += massI*molI->v();
                }
            }
        }
    }

    if(Pstream::parRun())
    {
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << mass << momentum;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalarField massProc;
                vectorField momProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> massProc >> momProc;
                }
    
                mass += massProc;
                momentum += momProc;
            }
        }
    }

    velocity_ = vector::zero;

    forAll(velocity_, n)
    {
        if(mass[n] > 0.0)
        {
            velocity_[n] = momentum[n]/mass[n];
        }
    }

}



void polyTemperatureBerendsenBinsNew::controlBeforeVelocityI()
{}

void polyTemperatureBerendsenBinsNew::controlBeforeMove()
{
}

void polyTemperatureBerendsenBinsNew::controlBeforeForces()
{}


void polyTemperatureBerendsenBinsNew::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyTemperatureBerendsenBinsNew::controlAfterForces()
{}

void polyTemperatureBerendsenBinsNew::controlAfterVelocityII()
{
    {
        Info << "polyTemperatureBerendsenBinsNew: measurement and control" << endl;

        const labelList& cells = mesh_.cellZones()[regionId_];

        const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];
                const vector& rI = molI->position();
    
                label n = binModel_->isPointWithinBin(rI, cell);
    
                if(n != -1)
                {
                    if(findIndex(molIds_, molI->id()) != -1)
                    {
//                         vector oldVel = molI->v();

                        if(!peculiar_)
                        {
                            if(componentControl_)
                            {
                                if(X_)
                                {
                                    molI->v().x() *= chi_[n];
                                }
                                if(Y_)
                                {
                                    molI->v().y() *= chi_[n];
                                }
                                if(Z_)
                                {
                                    molI->v().z() *= chi_[n];
                                }
                            }
                            else
                            {
                                molI->v() *= chi_[n];
                            }
                        }
                        else
                        {
                            if(angularControl_)
                            {
                                molI->v() = velocity_[n] + (chiLin_[n]*(molI->v() - velocity_[n]));
                                molI->pi() *= chiAng_[n];
                            }
                            else
                            {
                                molI->v() = velocity_[n] + (chi_[n]*(molI->v() - velocity_[n]));
                            }
                        }
                    }
                }
            }
        }
    }
}

void polyTemperatureBerendsenBinsNew::calculateProperties()
{

    const labelList& cells = mesh_.cellZones()[regionId_];

    {
        const List< DynamicList<polyMolecule*> >& cellOccupancy
            = molCloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];
                const vector& rI = molI->position();
    
                label n = binModel_->isPointWithinBin(rI, cell);
    
                if(n != -1)
                {
                    if((findIndex(molIds_, molI->id()) != -1) || measureFullTemperature_)
                    {
                        const scalar& massI = molCloud_.cP().mass(molI->id());
                        mass_[n] += massI;
                        momentum_[n] += massI*molI->v();
                    }
                }
            }
        }
    }

    {
        scalarField mass = mass_;
        vectorField momentum = momentum_;

        //- parallel communication
        if(Pstream::parRun())
        {
            //- sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << mass << momentum;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField massProc;
                    vectorField momProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> massProc >> momProc;
                    }
        
                    mass += massProc;
                    momentum += momProc;
                }
            }
        }

        velocity_ = vector::zero;

        forAll(velocity_, n)
        {
            if(mass[n] > 0.0)
            {
                velocity_[n] = momentum[n]/mass[n];
            }
        }

        //- reset 
        if(resetFieldsAtOutput_)
        {
            mass_ = 0.0;
            momentum_ = vector::zero;
        }
    }

    {
        const List< DynamicList<polyMolecule*> >& cellOccupancy
                                        = molCloud_.cellOccupancy();
    
        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<polyMolecule*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];
                const vector& rI = molI->position();
    
                label n = binModel_->isPointWithinBin(rI, cell);
    
                if(n != -1)
                {
                    if((findIndex(molIds_, molI->id()) != -1) || measureFullTemperature_)
                    {
                        const scalar& massI = molCloud_.cP().mass(molI->id());
    
                        kE_[n] += 0.5*massI*magSqr(molI->v() - velocity_[n]);
                        dof_[n] += molCloud_.cP().degreesOfFreedom(molI->id());
                        
                        const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));
                        const vector& molOmega(inv(molMoI) & molI->pi());
                        angularKe_[n] += 0.5*(molOmega & molMoI & molOmega);
                    }
                }
            }
        }
    }

//     if(time_.averagingTime())
    {
        scalarField kE = kE_;
        scalarField dof = dof_;
        scalarField angularKe = angularKe_;
        
       //- parallel processing
        if(Pstream::parRun())
        {
            //- sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << kE << angularKe << dof;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField kEProc;
                    scalarField angularKeProc;
                    scalarField dofProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> kEProc >>angularKeProc >> dofProc;
                    }

                    kE += kEProc;
                    angularKe += angularKeProc;
                    dof += dofProc;
                }
            }
        }

        measuredT_ = scalar(0.0);

        const scalar& kB = molCloud_.redUnits().kB();
         const scalar deltaTMD = time_.deltaT().value(); 
        chiAng_ = 1.0;
        chiLin_ = 1.0;
        
        forAll(measuredT_, n)
        {
            if(dof[n] > 0)
            {
                measuredT_[n] = (2.0*(kE[n]+ angularKe[n]))/(kB*dof[n]);

                chi_[n] = sqrt(1.0 + (deltaTMD/tauT_)*((temperature_/measuredT_[n]) - 1.0) );

                scalar measuredTLin = (2.0*kE[n])/(kB*dof[n]);
                scalar measuredTAng = (2.0*angularKe[n])/(kB*dof[n]);
                
                if(mag(measuredTLin) > 0)
                {
                    chiLin_[n] = sqrt(1.0 + (deltaTMD/tauT_)*((0.5*temperature_/measuredTLin) - 1.0) );
                }
                
                if(mag(measuredTAng) > 0)
                {
                    chiAng_[n] = sqrt(1.0 + (deltaTMD/tauT_)*((0.5*temperature_/measuredTAng) - 1.0) );
                }

                if(output_)
                {
                    Info << "measured T (R.U.) = " << measuredT_[n] 
                        << ", chi: " << chi_[n] 
                        << ", chiLin : " << chiLin_[n]
                        << ", chiAng : " << chiAng_[n]
                        << endl;
                }
            }
        }

         //- reset
        if(resetFieldsAtOutput_)
        {
            kE_ = 0.0;
            angularKe_ = 0.0;            
            dof_ = 0.0;
        }
    }
}



void polyTemperatureBerendsenBinsNew::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{

}



void polyTemperatureBerendsenBinsNew::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
/*
    if(propsDict_.found("tauT"))
    {
        tauT_ = readScalar(propsDict_.lookup("tauT"));
    }

    if (readStateFromFile_)
    {
        temperature_ = readScalar(propsDict_.lookup("temperature"));
    }*/

}




} // End namespace Foam

// ************************************************************************* //
