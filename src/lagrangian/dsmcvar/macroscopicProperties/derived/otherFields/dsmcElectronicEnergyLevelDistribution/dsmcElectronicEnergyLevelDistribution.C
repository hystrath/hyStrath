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

Measures the vibrational quantum level of every particle of the selected species in the selected zone.
The results are output in the time directory/uniform and are displayed with quantum level on the x-axis
and probability on the y-axis.

\*---------------------------------------------------------------------------*/

#include "dsmcElectronicEnergyLevelDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcElectronicEnergyLevelDistribution, 0);

addToRunTimeSelectionTable(dsmcField, dsmcElectronicEnergyLevelDistribution, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcElectronicEnergyLevelDistribution::dsmcElectronicEnergyLevelDistribution
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
    typeId_(-1),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    binWidth_(1), 
    distr_(binWidth_)
{
    word typeIdName = propsDict_.lookup("typeId");
    typeId_ = findIndex(cloud_.typeIdList(), typeIdName);

    if(typeId_ == -1)
    {
        FatalErrorIn("dsmcElectronicEnergyLevelDistribution::dsmcElectronicEnergyLevelDistribution()")
            << "Cannot find typeid: " << typeIdName << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcElectronicEnergyLevelDistribution::dsmcElectronicEnergyLevelDistribution()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcElectronicEnergyLevelDistribution::~dsmcElectronicEnergyLevelDistribution()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void dsmcElectronicEnergyLevelDistribution::createField()
{
    Info << "Initialising dsmcElectronicEnergyLevelDistribution field" << endl;
}


void dsmcElectronicEnergyLevelDistribution::calculateField()
{
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

                if(p->typeId() == typeId_)
                {
                    distr_.add(p->ELevel());
                }
            }
        }
    }
}

//- write field
void dsmcElectronicEnergyLevelDistribution::writeField()
{
    const Time& runTime = time_.time();

    if((runTime.outputTime()) && (time_.averagingTime()))
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List< Pair<scalar> > normalisedDistribution = distr_.normalised();
    
        if(Pstream::master())
        {
            scalarField xAxis (normalisedDistribution.size(), 0.0);
            scalarField yAxis (normalisedDistribution.size(), 0.0);

            forAll(normalisedDistribution, i)
            {
                xAxis[i] = normalisedDistribution[i].first() - (0.5*binWidth_);
                yAxis[i] = normalisedDistribution[i].second();
            }

            writeTimeData(timePath, "dsmcElectronicEnergyLevelDistribution_"+fieldName_+"_"+regionName_, xAxis, yAxis);

            if(time_.resetFieldsAtOutput())
            {
                distr_.clear();
            }
        }
    }
}

void dsmcElectronicEnergyLevelDistribution::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
