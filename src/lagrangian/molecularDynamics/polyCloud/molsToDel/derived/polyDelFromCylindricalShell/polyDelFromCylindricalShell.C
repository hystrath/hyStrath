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

#include "polyDelFromCylindricalShell.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromCylindricalShell, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromCylindricalShell, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromCylindricalShell::polyDelFromCylindricalShell
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    rOut_(readScalar(propsDict_.lookup("rOut"))),
    rIn_(readScalar(propsDict_.lookup("rIn")))
{

    // check if start point is in the mesh
   
    if(mesh_.findCell(startPoint_) == -1)
    {
        Info<< "WARNING: starting point " << startPoint_ 
            << " is selected outside the mesh."
            << endl;
    }

    if(mesh_.findCell(endPoint_) == -1)
    {
        Info<< "WARNING: end point " << endPoint_ 
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

polyDelFromCylindricalShell::~polyDelFromCylindricalShell()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelFromCylindricalShell::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    label initialSize = molCloud_.size();

    scalar rSEMag = mag(endPoint_ - startPoint_);

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        const vector& rI = mol().position();
        vector rSI = rI - startPoint_;
        scalar centreLineDistance = rSI & unitVector_;

        //- step 1: test polyMolecule is between starting point and end point
        if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
        {
            vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

            //step 2: test polyMolecule is within the radial part of the spherical shell
            if
            (
                (mag(pointOnCentreLine-rI) <= rOut_) &&
                (mag(pointOnCentreLine-rI) >= rIn_) 
            )
            {
                label molId = mol().id();

                if(findIndex(molIds_, molId) != -1)
                {
                    polyMolecule* molI = &mol();
                    molsToDel.append(molI);
                }
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
