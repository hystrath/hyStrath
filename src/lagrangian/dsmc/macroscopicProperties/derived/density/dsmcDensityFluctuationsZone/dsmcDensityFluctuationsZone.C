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

#include "dsmcDensityFluctuationsZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDensityFluctuationsZone, 0);

addToRunTimeSelectionTable(dsmcField, dsmcDensityFluctuationsZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDensityFluctuationsZone::dsmcDensityFluctuationsZone
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
    mass_(),
    USum_(),
    nParticles_(),
    timeIndex_(0),
    nInstantSteps_(0),
    binDensity_(),
    binVelocity_(),
    binDensityField_(),
    binVelocityField_()
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
            FatalErrorIn("dsmcDensityFluctuationsZone::dsmcDensityFluctuationsZone()")
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
        FatalErrorIn("dsmcDensityFluctuationsZone::dsmcDensityFluctuationsZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
   // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );
    
    scalar writeInterval = readScalar(t.controlDict().lookup("writeInterval"));
    scalar deltaT = t.deltaT().value();
    nInstantSteps_ = label(writeInterval/deltaT);

    const label& nBins = binModel_->nBins();

    if(Pstream::master())
    {
        mass_.setSize(nBins, 0.0);
        USum_.setSize(nBins, vector::zero);
        nParticles_.setSize(nBins, 0);
        binDensity_.setSize(nBins, 0.0);
        binVelocity_.setSize(nBins, 0.0);
        binDensityField_.setSize(nInstantSteps_);
        binVelocityField_.setSize(nInstantSteps_);
    }
    
//     timeIndex_ = 0;
    
//     scalar binWidth = binModel_->binPositions()[1]-binModel_->binPositions()[0];
//     domainLength_ = binWidth*nBins;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDensityFluctuationsZone::~dsmcDensityFluctuationsZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void dsmcDensityFluctuationsZone::createField()
{
    Info << "Initialising dsmcDensityFluctuationsZone field" << endl;
}


void dsmcDensityFluctuationsZone::calculateField()
{  
//     resetCounter_++;
   
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
                        
                        nParticles_[n] += cloud_.nParticle()*RWF;
                        USum_[n] += p->U();
                    }
                    else
                    {
                        nParticles_[n] += cloud_.nParticle();
                        USum_[n] += p->U();
                    }
                }
            }
        }
    }
    
    forAll(mass_, n)
    {
        scalar volume = binModel_->binVolume(n);
        
        binDensity_[n] = (nParticles_[n])/volume;
        binVelocity_[n] = USum_[n].x()/nParticles_[n];
    }

    binDensityField_[timeIndex_] = binDensity_;
    binVelocityField_[timeIndex_] = binVelocity_;

    timeIndex_++;
    
    mass_ = 0.0;
    nParticles_ = 0;
    USum_ = vector::zero;
    
//     const Time& runTime = time_.time();
    
//     if(runTime.outputTime())
//     {
//         resetCounter_ = 0;
//     }
}

//- write field
void dsmcDensityFluctuationsZone::writeField()
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
                "densityFluctuations_"+fieldName_+"_"+regionName_+".xy",
                timeField,
                binDensityField_,
                true
            );
            
            writeTimeData
            (
                casePath_,
                "velocityFluctuations_"+fieldName_+"_"+regionName_+".xy",
                timeField,
                binVelocityField_,
                true
            );
        }
    }
}

void dsmcDensityFluctuationsZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}

void dsmcDensityFluctuationsZone::setProperties()
{

}

} // End namespace Foam

// ************************************************************************* //
