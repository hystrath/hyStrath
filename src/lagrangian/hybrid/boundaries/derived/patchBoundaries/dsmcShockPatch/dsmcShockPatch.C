/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "dsmcShockPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcShockPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcShockPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcShockPatch::dsmcShockPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcShockPatch::~dsmcShockPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcShockPatch::initialConfiguration()
{}

void dsmcShockPatch::calculateProperties()
{}

void dsmcShockPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    if (p.classification() == 1000000)
    {
        td.keepParticle = false;
    }
    else
    {
        p.U().x() *= -1.0;
    }
}

void dsmcShockPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcShockPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}



} // End namespace Foam

// ************************************************************************* //
