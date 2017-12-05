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

#include "srbsFluctuationsComplexZone.H"
#include "addToRunTimeSelectionTable.H"
#include "complex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(srbsFluctuationsComplexZone, 0);

addToRunTimeSelectionTable(dsmcField, srbsFluctuationsComplexZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
srbsFluctuationsComplexZone::srbsFluctuationsComplexZone
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    resetCounter_(0.0),
    domainLength_(),
    Rt_(),
    RtStar_(),
    RtField_(time_.totalNAvSteps()+1),
    RtStarField_(time_.totalNAvSteps()+1)
{

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
            FatalErrorIn("srbsFluctuationsComplexZone::srbsFluctuationsComplexZone()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------


    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("srbsFluctuationsComplexZone::srbsFluctuationsComplexZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
   // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );

    const label& nBins = binModel_->nBins();
    n_.setSize(nBins, 0.0);
    rtSin_.setSize(nBins, 0.0);
    rtCos_.setSize(nBins, 0.0);

    scalar writeInterval = readScalar(t.controlDict().lookup("writeInterval"));
    scalar deltaT = t.deltaT().value();
    nInstantSteps_ = label(writeInterval/deltaT);

    if(Pstream::master())
    {
        RtField_.setSize(nInstantSteps_);
        RtStarField_.setSize(nInstantSteps_);
    }
    
    timeIndex_ = 0;
    
    scalar binWidth = binModel_->binPositions()[1]-binModel_->binPositions()[0];
    domainLength_ = binWidth*nBins;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

srbsFluctuationsComplexZone::~srbsFluctuationsComplexZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void srbsFluctuationsComplexZone::createField()
{
    Info << "Initialising srbsFluctuationsComplexZone field" << endl;
}


void srbsFluctuationsComplexZone::calculateField()
{
       
    const List< DynamicList<dsmcParcel*> >& cellOccupancy
        = cloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];
    
    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cellI];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            
            const vector& rI = p->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    if(cloud_.axisymmetric())
                    {
                        const point& cC = cloud_.mesh().cellCentres()[cellI];
                        scalar radius = cC.y();
                        
                        scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                        
                        n_[n] += cloud_.nParticle()*RWF;
                    }
                    else
                    {
                        n_[n] += cloud_.nParticle();
                    }
                }
            }
        }
    }

    Rt_.Re() = 0.0;
    Rt_.Im() = 0.0;
    
    RtStar_.Re() = 0.0;
    RtStar_.Im() = 0.0;
    
    const scalarField& binCentres = binModel_->binPositions();
    
    forAll(binCentres, n)
    {
        scalar volume = binModel_->binVolume(n);
        rtSin_[n] = (n_[n]/volume)*sin(twoPi*binCentres[n]/domainLength_);
        rtCos_[n] = (n_[n]/volume)*cos(twoPi*binCentres[n]/domainLength_);
        
        Rt_.Re() += rtCos_[n];
        Rt_.Im() += rtSin_[n];
        
        RtStar_.Re() += rtCos_[n];
        RtStar_.Im() -= rtSin_[n];
    }
    
    Rt_ /= binModel_->nBins();
    RtStar_ /= binModel_->nBins();

    RtField_[timeIndex_] = Rt_;
    RtStarField_[timeIndex_] = RtStar_;
    timeIndex_++;
    n_ = scalar(0.0);
}


//- write field
void srbsFluctuationsComplexZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        timeIndex_ = 0;
        
        if(Pstream::master())
        {
            fileName timePath(runTime.path()/runTime.timeName()/"uniform");
        
            if (!isDir(timePath))
            {
                mkDir(timePath);
            }
            
            scalar dt = time_.time().deltaT().value();
            scalarField timeField (nInstantSteps_, 0.0);

            forAll(timeField, t)
            {
                timeField[nInstantSteps_-t-1] = time_.time().timeOutputValue()-dt*t;
            }
            
            writeTimeData
            (
                casePath_,
                "srbsFluctuations_"+fieldName_+"_"+regionName_+"_Rt.xy",
                timeField,
                RtField_,
                true
            );
            
            writeTimeData
            (
                casePath_,
                "srbsFluctuations_"+fieldName_+"_"+regionName_+"_RtStar.xy",
                timeField,
                RtStarField_,
                true
            );
        }
    }
}



void srbsFluctuationsComplexZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}

void srbsFluctuationsComplexZone::setProperties()
{

}

} // End namespace Foam

// ************************************************************************* //
