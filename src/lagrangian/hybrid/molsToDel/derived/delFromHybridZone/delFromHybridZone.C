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

#include "delFromHybridZone.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(delFromHybridZone, 0);

addToRunTimeSelectionTable
(
    molsToDeleteModel,
    delFromHybridZone,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void delFromHybridZone::findMolsToDel()
{
    const labelList& cells = mesh_.cellZones()[regionId_];

    DynamicList<dsmcParcel*> molsToDel;

    label initialSize = cloud_.size();

    /*forAllIter(dsmcCloud, cloud_, p)
    {
        const label& cellI =  p().cell();

        if(findIndex(cells, cellI) != -1)
        {
            if(findIndex(typeIds_, p().typeId()) != -1)
            {
                dsmcParcel* pI = &p();
                molsToDel.append(pI);
            }
        }
    }*/

    forAll(cells, cellsI)
    {
        const label cell = cells[cellsI];
        molsToDel.append(cloud_.cellOccupancy()[cell]);
    }

    molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    const label molsKept = initialSize - molsToDel.size();

    Info<< "  Zone: " << regionName_ << tab
        << "particles removed: " << molsToDel.size()
        << ", new cloud size: " << molsKept
        << endl;


    // as a precaution: rebuild cell occupancy
//     cloud_.rebuildCellOccupancy();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
delFromHybridZone::delFromHybridZone
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    molsToDeleteModel(cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    typeIds_()
{
    Info<< propsDict_.name() << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("delFromHybridZone::delFromHybridZone(...)")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << mesh_.time().system()/"molsToDeleteDict"
            << exit(FatalError);
    }


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
            FatalErrorIn("delFromHybridZone::delFromHybridZone(...)")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"molsToDeleteDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

delFromHybridZone::~delFromHybridZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void delFromHybridZone::update()
{
    findMolsToDel();
}


} // End namespace Foam

// ************************************************************************* //
