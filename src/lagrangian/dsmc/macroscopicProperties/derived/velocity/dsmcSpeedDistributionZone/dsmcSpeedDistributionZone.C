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

#include "dsmcSpeedDistributionZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcSpeedDistributionZone, 0);

addToRunTimeSelectionTable(dsmcField, dsmcSpeedDistributionZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcSpeedDistributionZone::dsmcSpeedDistributionZone
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
    typeIds_(),
    timeDictVel_(dict.subDict("timePropertiesForVelocity")),
    timeVel_(t, timeDictVel_),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    UMean_(vector::zero),
    Ucollected_(vector::zero),
    nParcels_(0),
    binWidth_(readScalar(propsDict_.lookup("binWidth"))), 
    distr_(binWidth_)
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
            FatalErrorIn("dsmcSpeedDistributionZone::dsmcSpeedDistributionZone()")
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
        FatalErrorIn("dsmcSpeedDistributionZone::dsmcSpeedDistributionZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSpeedDistributionZone::~dsmcSpeedDistributionZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void dsmcSpeedDistributionZone::createField()
{
    Info << "Initialising dsmcSpeedDistributionZone field" << endl;
}


void dsmcSpeedDistributionZone::calculateField()
{
    timeVel_++;

    if(timeVel_.samplingTime())
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

                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    nParcels_++;
                    Ucollected_ += p->U();
                }
            }
        }
    }

    if(timeVel_.averagingTime())
    {
//         resetIndex_++;

        vector Ucollected = Ucollected_;
        label nParcels = nParcels_;

        if(Pstream::parRun())
        {
            reduce(Ucollected, sumOp<vector>());
            reduce(nParcels, sumOp<label>());
        }

        if(nParcels_ > 0)
        {
            UMean_ = Ucollected/nParcels;
            Info << "dsmcSpeedDistributionZone, averaged velocity " << UMean_ << endl;
        }

        if(time_.resetFieldsAtOutput())
        {
            Ucollected_ = vector::zero;
            nParcels_ = 0;
        }
    }

    if(time_.samplingTime())
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

                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    distr_.add(mag(p->U() - UMean_));

//                     Info << " adding: "  << mag(p->U() - UMean_) << endl;
                }
            }
        }
    }
}


//- return field
// const volScalarField& dsmcSpeedDistributionZone::densityField() const
// {
//     return rhoN_;
// }

//- write field
void dsmcSpeedDistributionZone::writeField()
{
    const Time& runTime = time_.time();

    if((runTime.outputTime()) && (time_.averagingTime()))
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List< Pair<scalar> > normalisedDistriubtion = distr_.normalised();

        label nSize = normalisedDistriubtion.size();

        if (Pstream::parRun())
        {
            reduce(nSize, sumOp<label>());
        }

        scalarField xAxis (nSize, 0.0);
        scalarField yAxis (nSize, 0.0);

        forAll(normalisedDistriubtion, i)
        {
            xAxis[i] = normalisedDistriubtion[i].first();
            yAxis[i] = normalisedDistriubtion[i].second();
        }

        if (Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << xAxis << yAxis;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField xAxisProc;
                    scalarField yAxisProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> xAxisProc >> yAxisProc;
                    }
    
                    forAll(xAxisProc, i)
                    {
                        xAxis[i] += xAxisProc[i];
                        yAxis[i] += yAxisProc[i];
                    }
                }
            }
        }

        writeTimeData(timePath, "speedDistribution_"+fieldName_+"_"+regionName_, xAxis, yAxis);
        
        if(time_.resetFieldsAtOutput())
        {
            distr_.clear();
        }
    }
}

void dsmcSpeedDistributionZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
}

// const propertyField& dsmcSpeedDistributionZone::fields() const
// {
//     return  fields_;
// }




} // End namespace Foam

// ************************************************************************* //
