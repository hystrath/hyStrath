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

Controls translational temperature on a zone basis in a user-defined zone in a mesh.
Allows for better statistics than a cell-by-cell, but resolution is lost.

\*---------------------------------------------------------------------------*/

#include "temperatureZoneController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(temperatureZoneController, 0);

addToRunTimeSelectionTable(dsmcStateController, temperatureZoneController, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temperatureZoneController::temperatureZoneController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    timeDictVel_(propsDict_.subDict("timePropertiesForVelocity")),
    timeVel_(t, timeDictVel_),
    tauT_(readScalar(propsDict_.lookup("tauT"))),
    componentControl_(false),
    X_(false),
    Y_(false),
    Z_(false),
    typeIds_(),
    massV_(scalar(0.0)),
    momV_(vector::zero),
    UMean_(vector::zero),
    mcc_(scalar(0.0)),
    m_(scalar(0.0)),
    nParcels_(scalar(0.0)),
    measuredTranslationalTemperature_(scalar(0.0)),
    chi_(scalar(0.0))
{
    if (propsDict_.found("componentControl"))
    {
        componentControl_ = Switch(propsDict_.lookup("componentControl"));

        if(componentControl_)
        {
            if (propsDict_.found("X"))
            {
                X_ = true;
            }
    
            if (propsDict_.found("Y"))
            {
                Y_ = true;
            }
    
            if (propsDict_.found("Z"))
            {
                Z_ = true;
            }
            
            Info << "X_ = " << X_ << ", Y_ = " << Y_ << ", Z = " << Z_ << endl;

            if(!X_ && !Y_ && !Z_)
            {
                FatalErrorIn("temperatureControllerBins::temperatureControllerBins()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }


    setProperties();

    measuredTranslationalTemperature_ = temperature_;
    
    // standard to reading typeIds ------------ 
    const List<word> molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcInflowPatch::dsmcInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureZoneController::~temperatureZoneController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureZoneController::initialConfiguration()
{}

void temperatureZoneController::calculateProperties()
{
    timeVel_++;

    const labelList& cells = mesh_.cellZones()[regionId_];

    if(timeVel_.samplingTime())
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
    
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    scalar nParticle = cloud_.nParticle();
                    
                    const scalar& RWF = cloud_.getRWF_cell(cellI);
                    
                    nParticle *= RWF;
                    
                    const scalar mass = cloud_.constProps(p->typeId()).mass()*nParticle;
                    
                    massV_ += mass;
                    momV_ += p->U()*mass;
                }
            }
        }
    }

    if(timeVel_.averagingTime())
    {
        scalar massV = massV_;
        vector momV = momV_;

        //- parallel communication
        if(Pstream::parRun())
        {
            // sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << massV << momV;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar massVProc;
                    vector momVProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> massVProc >> momVProc;
                    }
    
                    massV += massVProc;
                    momV += momVProc;
                }
            }
        }
            
        UMean_ = vector::zero;

        if(massV > 0.0)
        {
        UMean_ = momV/massV;
        }

        //- reset 
        if(time_.resetFieldsAtOutput())
        {
            massV_ = 0.0;
            momV_ = vector::zero;
        }
    }
    
    if(time_.samplingTime())
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
    
        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    scalar nParticle = cloud_.nParticle();
                    
                    const scalar& RWF = cloud_.getRWF_cell(cell);
                    
                    nParticle *= RWF;
                    
                    const scalar mass = cloud_.constProps(p->typeId()).mass()*nParticle;
                    
                    mcc_ += mass*mag(p->U())*mag(p->U());
                    m_ += mass;
                    nParcels_ += nParticle;
                }
            }
        }
    }

    if(time_.averagingTime())
    {
        scalar mcc = mcc_;
        scalar m = m_;
        scalar nParcels = nParcels_;

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
                        toNeighbour << mcc << m << nParcels;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar mccProc;
                    scalar mProc;
                    scalar nParcelsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> mccProc >> mProc >> nParcelsProc;
                    }

                    mcc += mccProc;
                    m += mProc;
                    nParcels += nParcelsProc;
                }
            }
        }
    
        measuredTranslationalTemperature_ = scalar(0.0);
        
        const scalar& deltaTDSMC = mesh_.time().deltaTValue(); // time step 

        if(nParcels > 0.0)
        {
            measuredTranslationalTemperature_ = (1.0/(3.0*physicoChemical::k.value()))
                *(
                        ((mcc/(nParcels)))
                        - ((m/(nParcels))*mag(UMean_)*mag(UMean_))
                );

            chi_ = sqrt(1.0 + (deltaTDSMC/tauT_)*((temperature_/measuredTranslationalTemperature_) - 1.0) );
            
            Info<< "target temperature: " << temperature_ 
                    << " UMean_ : " << UMean_
                    << " measured T: " << measuredTranslationalTemperature_
                    << " chi: " << chi_
                    << endl;
        }

         //- reset
        if(time_.resetFieldsAtOutput())
        {
            mcc_ = 0.0;
            m_ = 0.0;
            nParcels_ = 0.0;
        }
    }
}

void temperatureZoneController::controlParcelsBeforeMove()
{
    if(control_ && time_.controlTime())
    {
        Info << "temperatureController: control" << endl;
        
        const labelList& cells = mesh_.cellZones()[regionId_];

        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    if(componentControl_)
                    {
                        if(X_ && chi_ > 0)
                        {
                                p->U().x() *= chi_;
                        }
                        if(Y_ && chi_ > 0)
                        {
                                p->U().y() *= chi_;
                        }
                        if(Z_ && chi_ > 0)
                        {
                                p->U().z() *= chi_;
                        }
                    }
                    else
                    {
                        if(chi_ > 0)
                        {
                                p->U() -= UMean_;
                                p->U() *= chi_;
                                p->U() += UMean_;
                        }
                    }
                }
            }
        }
    }
}


void temperatureZoneController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
}


void temperatureZoneController::controlParcelsBeforeCollisions()
{}

void temperatureZoneController::controlParcelsAfterCollisions()
{}

void temperatureZoneController::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

void temperatureZoneController::setProperties()
{
    temperature_ = readScalar(propsDict_.lookup("controlTemperature"));
}


}
// End namespace Foam

// ************************************************************************* //
