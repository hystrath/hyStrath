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

Controls translational temperature in bins in a user-defined zone in a mesh.
Allows for a mixture of cell-by-cell and zone treatment, improved statistics
compared to cell-by-cell and increased resolution compared to zone treatment.

\*---------------------------------------------------------------------------*/

#include "temperatureBinsController.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(temperatureBinsController, 0);

addToRunTimeSelectionTable(dsmcStateController, temperatureBinsController, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temperatureBinsController::temperatureBinsController
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
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    nBins_(readLabel(propsDict_.lookup("nBins"))),
    binWidth_(mag(endPoint_ - startPoint_)/(nBins_)),
    startTemperature_(readScalar(propsDict_.lookup("startTemperature"))),
    endTemperature_(readScalar(propsDict_.lookup("endTemperature"))),
    binTemperatures_(nBins_, 0.0),
    typeIds_(),
    magRadii_(nBins_, 0.0),
    massV_(nBins_, scalar(0.0)),
    momV_(nBins_, vector::zero),
    UMean_(nBins_, vector::zero),
    mcc_(nBins_, scalar(0.0)),
    m_(nBins_, scalar(0.0)),
    nParcels_(nBins_, scalar(0.0)),
    measuredTranslationalTemperature_(nBins_, scalar(0.0)),
    chi_(nBins_, 0.0)
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

            if(!X_ && !Y_ && !Z_)
            {
                FatalErrorIn("temperatureControllerBins::temperatureControllerBins()")
                    << "At least one component (X, Y, Z) should be chosen " << nl << "in: "
                    << time_.time().system()/"controllersDict"
                    << exit(FatalError);
            }
        }
    }

//     readProperties();
	
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

temperatureBinsController::~temperatureBinsController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureBinsController::initialConfiguration()
{
    const scalar& deltaTBin = (endTemperature_ - startTemperature_)/(nBins_ - 1);
    
    forAll(binTemperatures_, n)
    {
        binTemperatures_[n] = startTemperature_ + n*deltaTBin;
    }
    
    Info << "binTemperatures_ = " << binTemperatures_ << endl;
}

void temperatureBinsController::calculateProperties()
{
    timeVel_++;

    scalar rSEMag = mag(endPoint_ - startPoint_);

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

                const vector& rI = p->position();
                vector rSI = rI - startPoint_;
                scalar rD = rSI & unitVector_;
                label n = label(rD/binWidth_);

                if((n == nBins_) && (rD <= rSEMag))
                {
                        n--;
                }

                if
                (
                    (rD >= 0) &&
                    (n >= 0) &&
                    (n < nBins_)
                )
                {
                    if(findIndex(typeIds_, p->typeId()) != -1)
                    {

                        const scalar nParticle = cloud_.nParticles(cellI);
                    
                        const scalar mass = cloud_.constProps(p->typeId()).mass()*nParticle;
                        
                        massV_[n] += mass;
                        momV_[n] += p->U()*mass;
                    }
                }
            }
        }
    }

    if(timeVel_.averagingTime())
    {
        scalarField mass = massV_;
        vectorField mom = momV_;

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
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << mass << mom;
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
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> massProc >> momProc;
                    }
    
                    mass += massProc;
                    mom += momProc;
                }
            }
        }

        UMean_ = vector::zero;

        forAll(UMean_, n)
        {
            if(mass[n] > 0.0)
            {
                UMean_[n] = mom[n]/mass[n];
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
    
        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];

                const vector& rI = p->position();
                vector rSI = rI - startPoint_;
                scalar rD = rSI & unitVector_;
                label n = label(rD/binWidth_);

                if((n == nBins_) && (rD <= rSEMag))
                {
                    n--;
                }

                if
                (
                    (rD >= 0) &&
                    (n >= 0) &&
                    (n < nBins_)
                )
                {
                    if(findIndex(typeIds_, p->typeId()) != -1)
                    {
                        const scalar nParticle = cloud_.nParticles(cell);
                    
                        const scalar mass = cloud_.constProps(p->typeId()).mass()*nParticle;
                        
                        mcc_[n] += mass*mag(p->U())*mag(p->U());
                        m_[n] += mass;
                        nParcels_[n] += nParticle;
                    }
                }
            }
        }
    }

    if(time_.averagingTime())
    {
        scalarField mcc = mcc_;
        scalarField m = m_;
        scalarField nParcels = nParcels_;

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
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << mcc << m << nParcels;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField mccProc;
                    scalarField mProc;
                    scalarField nParcelsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
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

        forAll(measuredTranslationalTemperature_, n)
        {
            if(nParcels[n] > 0.0)
            {
                measuredTranslationalTemperature_[n] = (1.0/(3.0*physicoChemical::k.value()))
                    *(
                            ((mcc[n]/(nParcels[n])))
                            - ((m[n]/(nParcels[n]))*mag(UMean_[n])*mag(UMean_[n]))
                    );

                chi_[n] = sqrt(1.0 + (deltaTDSMC/tauT_)*((binTemperatures_[n]/measuredTranslationalTemperature_[n]) - 1.0) );
                
                Info<< "target temperature: " << binTemperatures_[n] 
                    << " measured T: " << measuredTranslationalTemperature_[n]
                    << " chi: " << chi_[n]
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

void temperatureBinsController::controlParcelsBeforeMove()
{
	
}


void temperatureBinsController::output
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


void temperatureBinsController::controlParcelsBeforeCollisions()
{
    if(control_ && time_.controlTime())
    {
        Info << "temperatureBinsController: control" << endl;
        
        const labelList& cells = mesh_.cellZones()[regionId_];

        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        scalar rSEMag = mag(endPoint_ - startPoint_);

        forAll(cells, c)
        {
            const label& cell = cells[c];
            const List<dsmcParcel*>& molsInCell = cellOccupancy[cell];
    
            forAll(molsInCell, mIC)
            {
                dsmcParcel* p = molsInCell[mIC];
                
                const vector& rI = p->position(); //position
                vector rSI = rI - startPoint_;
                scalar rD = rSI & unitVector_;
                label n = label(rD/binWidth_);

                if((n == nBins_) && rD <= rSEMag)
                {
                    n--;
                }

                if
                (
                    (rD >= 0) &&
                    (n >= 0) &&
                    (n < magRadii_.size())
                )
                {
                    if(findIndex(typeIds_, p->typeId()) != -1)
                    {
                        if(componentControl_)
                        {
                            if(X_ && chi_[n] > 0)
                            {
                                p->U().x() *= chi_[n];
                            }
                            if(Y_ && chi_[n] > 0)
                            {
                                p->U().y() *= chi_[n];
                            }
                            if(Z_ && chi_[n] > 0)
                            {
                                p->U().z() *= chi_[n];
                            }
                        }
                        else
                        {
                            if(chi_[n] > 0)
                            {
                                p->U() -= UMean_[n];
                                p->U() *= chi_[n];
                                p->U() += UMean_[n];
                            }
                        }
                    }
                }
            }
        }
    }
}

void temperatureBinsController::controlParcelsAfterCollisions()
{}

void temperatureBinsController::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

}

void temperatureBinsController::setProperties()
{

}


}
// End namespace Foam

// ************************************************************************* //
