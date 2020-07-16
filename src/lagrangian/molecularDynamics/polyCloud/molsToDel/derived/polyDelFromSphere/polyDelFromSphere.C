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

#include "polyDelFromSphere.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromSphere, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromSphere, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromSphere::polyDelFromSphere
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    radius_(readScalar(propsDict_.lookup("radius")))
{

    // check if start point is in the mesh

    if(mesh_.findCell(startPoint_) == -1)
    {
        Info<< "WARNING: starting point " << startPoint_
            << " is selected outside the mesh."
            << endl;
    }


    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelFromSphere::~polyDelFromSphere()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelFromSphere::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    label initialSize = molCloud_.size();

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        const vector& rI = mol().position();

        scalar distance = mag(rI - startPoint_);

        if(distance <= radius_)
        {
            label molId = mol().id();

            if(findIndex(molIds_, molId) != -1)
            {
                polyMolecule* molI = &mol();
                molsToDel.append(molI);
            }
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
