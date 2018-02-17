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

#include "dsmcShockGeneralPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcShockGeneralPatch, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcShockGeneralPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcShockGeneralPatch::dsmcShockGeneralPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    boundaryXPosition_(),
    boundaryXVelocity_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcShockGeneralPatch::~dsmcShockGeneralPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcShockGeneralPatch::initialConfiguration()
{}

void dsmcShockGeneralPatch::calculateProperties()
{

}

void dsmcShockGeneralPatch::controlParcelsBeforeMove()
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    const scalar& boundaryXPosition = boundaryXPosition_;
    const scalar& boundaryXVelocity = boundaryXVelocity_;

    const List<DynamicList<dsmcParcel*> >& cellOccupancy
        = cloud_.cellOccupancy();

    forAll(cells_, cellI)
    {
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[cellI]];
                    
        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            if(p->U().x() * deltaT
                > 2.0 * boundaryXVelocity * deltaT + boundaryXPosition
                - p->position().x())
            {
                p->U().x() -= 2.0 * boundaryXVelocity;
            }
            else if(p->position().x() + p->U().x() * deltaT
                > boundaryXPosition)
            {
                // classification should match the one in dsmcShockPatch.C
                p->classification() = 1000000;
            }
        }
    }
}

void dsmcShockGeneralPatch::controlParcelsBeforeCollisions()
{

}

void dsmcShockGeneralPatch::controlParcelsAfterCollisions()
{
}

void dsmcShockGeneralPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcShockGeneralPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

//     setProperties();
}



void dsmcShockGeneralPatch::setProperties()
{
    boundaryXPosition_ = readScalar(propsDict_.lookup("boundaryXPosition"));
    boundaryXVelocity_ = readScalar(propsDict_.lookup("boundaryXVelocity"));

    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcShockGeneralPatch::dsmcShockGeneralPatch()")
            << "Cannot have zero typeIds being inserd." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

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

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcShockGeneralPatch::dsmcShockGeneralPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
}



} // End namespace Foam

// ************************************************************************* //
