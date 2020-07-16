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

#include "polyDelFromMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromMesh, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromMesh, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromMesh::polyDelFromMesh
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelFromMesh::~polyDelFromMesh()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelFromMesh::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

//     label counter = 0;
    label initialSize = molCloud_.size();

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        label molId = mol().id();

        if(findIndex(molIds_, molId) != -1)
        {
            polyMolecule* molI = &mol();
            molsToDel.append(molI);
        }
    }

    //molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size()
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
}


} // End namespace Foam

// ************************************************************************* //
