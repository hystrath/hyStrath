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

#include "polySpecularPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polySpecularPatch, 0);

addToRunTimeSelectionTable(polyPatchBoundary, polySpecularPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polySpecularPatch::polySpecularPatch
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyPatchBoundary(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polySpecularPatch::~polySpecularPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polySpecularPatch::initialConfiguration()
{}

void polySpecularPatch::calculateProperties()
{}

void polySpecularPatch::controlMol
(
    polyMolecule& mol,
    polyMolecule::trackingData& td
)
{
//     Info << "molecule at pos = " << mol.position() << endl;

    const label& faceI = mol.face();

    vector nF = mesh_.faceAreas()[faceI];
    nF /= mag(nF);

    scalar Un = mol.v() & nF;

    if (Un > 0.0)
    {
        mol.v() -= 2.0*Un*nF;
    }
}

void polySpecularPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void polySpecularPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}



} // End namespace Foam

// ************************************************************************* //
