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

Controls translational temperature on a cell-by-cell basis in a user-defined zone on a mesh.
Provides the best resolution, but the worst statistics (low number of particles per cell).
No need for parallel processing as everything is on a per-cell basis.

\*---------------------------------------------------------------------------*/

#include "temperatureController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(temperatureController, 0);

addToRunTimeSelectionTable(dsmcStateController, temperatureController, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temperatureController::temperatureController
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
    massV_(controlZone().size(), scalar(0.0)),
    momV_(controlZone().size(), vector::zero),
    UMean_(controlZone().size(), vector::zero),
    mcc_(controlZone().size(), scalar(0.0)),
    m_(controlZone().size(), scalar(0.0)),
    nParcels_(controlZone().size(), scalar(0.0)),
    measuredTranslationalTemperature_(controlZone().size(), scalar(0.0)),
    chi_(controlZone().size(), 0.0)
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

temperatureController::~temperatureController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureController::initialConfiguration()
{
	
}

void temperatureController::calculateProperties()
{
    timeVel_++;

    const labelList& cells = mesh_.cellZones()[regionId_];

    if(timeVel_.samplingTime())
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        forAll(controlZone(), c)
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
                    
                    massV_[c] += mass;
                    momV_[c] += p->U()*mass;
                }
            }
        }
    }

    if(timeVel_.averagingTime())
    {
        UMean_ = vector::zero;

        forAll(UMean_, c)
        {
            if(massV_[c] > 0.0)
            {
                UMean_[c] = momV_[c]/massV_[c];
            }
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
    
        forAll(controlZone(), c)
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
                    
                    mcc_[c] += mass*mag(p->U())*mag(p->U());
                    m_[c] += mass;
                    nParcels_[c] += nParticle;
                }
            }
        }
    }

    if(time_.averagingTime())
    {
        measuredTranslationalTemperature_ = scalar(0.0);
            
            const scalar& deltaTDSMC = mesh_.time().deltaTValue(); // time step 

        forAll(measuredTranslationalTemperature_, c)
        {
            if(nParcels_[c] > 0.0)
            {
                measuredTranslationalTemperature_[c] = (1.0/(3.0*physicoChemical::k.value()))
                        *(
                                ((mcc_[c]/(nParcels_[c])))
                                - ((m_[c]/(nParcels_[c]))*mag(UMean_[c])*mag(UMean_[c]))
                        );

                chi_[c] = sqrt(1.0 + (deltaTDSMC/tauT_)*((temperature_/measuredTranslationalTemperature_[c]) - 1.0) );
                
                Info<< "target temperature: " << temperature_ 
                    << " UMean_ : " << UMean_[c] 
                    << " measured T: " << measuredTranslationalTemperature_[c]
                    << " chi: " << chi_[c]
                    << endl;
            }
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

void temperatureController::controlParcelsBeforeMove()
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
                        if(X_ && chi_[c] > 0)
                        {
                            p->U().x() *= chi_[c];
                        }
                        if(Y_ && chi_[c] > 0)
                        {
                            p->U().y() *= chi_[c];
                        }
                        if(Z_ && chi_[c] > 0)
                        {
                            p->U().z() *= chi_[c];
                        }
                    }
                    else
                    {
                        if(chi_[c] > 0)
                        {
                            p->U() -= UMean_[c];
                            p->U() *= chi_[c];
                            p->U() += UMean_[c];
                        }
                    }
                }
            }
        }
    }
}


void temperatureController::output
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


void temperatureController::controlParcelsBeforeCollisions()
{
    
}

void temperatureController::controlParcelsAfterCollisions()
{}

void temperatureController::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

}

void temperatureController::setProperties()
{
    temperature_ = readScalar(propsDict_.lookup("controlTemperature"));
}


}
// End namespace Foam

// ************************************************************************* //
