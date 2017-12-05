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

Outputs the velocity distribution function in x, y and z components in a user defined zone.

\*---------------------------------------------------------------------------*/

#include "dsmcVelocityDistributionZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVelocityDistributionZone, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVelocityDistributionZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVelocityDistributionZone::dsmcVelocityDistributionZone
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
    distrX_(binWidth_),
    distrY_(binWidth_),
    distrZ_(binWidth_)
    
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
            FatalErrorIn("dsmcVelocityDistributionZone::dsmcVelocityDistributionZone()")
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
        FatalErrorIn("dsmcVelocityDistributionZone::dsmcVelocityDistributionZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVelocityDistributionZone::~dsmcVelocityDistributionZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void dsmcVelocityDistributionZone::createField()
{
    Info << "Initialising dsmcVelocityDistributionZone field" << endl;
}


void dsmcVelocityDistributionZone::calculateField()
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
//                     if(cloud_.axisymmetric())
//                     {
//                         const point& cC = cloud_.mesh().cellCentres()[cellI];
//                         scalar radius = cC.y();
//                         
//                         scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
//                         
//                         nParcels_ += RWF;
//                         Ucollected_ += p->U();
//                     }
//                     else
//                     {
                        nParcels_++;
                        Ucollected_ += p->U();
//                     }
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
            Info << "dsmcVelocityDistributionZone, averaged velocity " << UMean_ << endl;
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
                    distrX_.add((p->U().x() - UMean_.x()));
                    distrY_.add((p->U().y() - UMean_.y()));
                    distrZ_.add((p->U().z() - UMean_.z()));
//                     Info << " adding: "  << mag(p->U() - UMean_) << endl;
                }
            }
        }
    }
}


//- return field
// const volScalarField& dsmcVelocityDistributionZone::densityField() const
// {
//     return rhoN_;
// }

//- write field
void dsmcVelocityDistributionZone::writeField()
{
    const Time& runTime = time_.time();

    if((runTime.outputTime()) && (time_.averagingTime()))
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List< Pair<scalar> > rawDistriubtionX = distrX_.raw();
        List< Pair<scalar> > rawDistriubtionY = distrY_.raw();
        List< Pair<scalar> > rawDistriubtionZ = distrZ_.raw();

        label nSizeX = rawDistriubtionX.size();
        label nSizeY = rawDistriubtionY.size();
        label nSizeZ = rawDistriubtionZ.size();

        if (Pstream::parRun())
        {
            reduce(nSizeX, sumOp<label>());
            reduce(nSizeY, sumOp<label>());
            reduce(nSizeZ, sumOp<label>());
        }

        scalarField xAxisX (nSizeX, 0.0);
        scalarField yAxisX (nSizeX, 0.0);
        scalarField xAxisY (nSizeY, 0.0);
        scalarField yAxisY (nSizeY, 0.0);
        scalarField xAxisZ (nSizeZ, 0.0);
        scalarField yAxisZ (nSizeZ, 0.0);
        

        forAll(rawDistriubtionX, i)
        {
            xAxisX[i] = rawDistriubtionX[i].first();
            yAxisX[i] = rawDistriubtionX[i].second();
        }
        
        forAll(rawDistriubtionY, i)
        {
            xAxisY[i] = rawDistriubtionY[i].first();
            yAxisY[i] = rawDistriubtionY[i].second();
        }
        
        forAll(rawDistriubtionZ, i)
        {
            xAxisZ[i] = rawDistriubtionZ[i].first();
            yAxisZ[i] = rawDistriubtionZ[i].second();
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
                        toNeighbour << xAxisX << yAxisX << xAxisY 
                                    << yAxisY << xAxisZ << yAxisZ;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField xAxisXProc;
                    scalarField yAxisXProc;
                    scalarField xAxisYProc;
                    scalarField yAxisYProc;
                    scalarField xAxisZProc;
                    scalarField yAxisZProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> xAxisXProc >> yAxisXProc >> xAxisYProc 
                                    >> yAxisYProc >> xAxisZProc >> yAxisZProc;
                    }
    
                    forAll(xAxisXProc, i)
                    {
                        xAxisX[i] += xAxisXProc[i];
                        yAxisX[i] += yAxisXProc[i];
                    }
                    
                    forAll(xAxisYProc, i)
                    {
                        xAxisY[i] += xAxisYProc[i];
                        yAxisY[i] += yAxisYProc[i];
                    }
                    
                    forAll(xAxisZProc, i)
                    {
                        xAxisZ[i] += xAxisZProc[i];
                        yAxisZ[i] += yAxisZProc[i];
                    }
                }
            }
        }

        writeTimeData(timePath, "velocityDistributionX_"+fieldName_+"_"+regionName_, xAxisX, yAxisX);
        writeTimeData(timePath, "velocityDistributionY_"+fieldName_+"_"+regionName_, xAxisY, yAxisY);
        writeTimeData(timePath, "velocityDistributionZ_"+fieldName_+"_"+regionName_, xAxisZ, yAxisZ);
        
        if(time_.resetFieldsAtOutput())
        {
            distrX_.clear();
            distrY_.clear();
            distrZ_.clear();
        }
    }
}

void dsmcVelocityDistributionZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
}

// const propertyField& dsmcVelocityDistributionZone::fields() const
// {
//     return  fields_;
// }




} // End namespace Foam

// ************************************************************************* //
