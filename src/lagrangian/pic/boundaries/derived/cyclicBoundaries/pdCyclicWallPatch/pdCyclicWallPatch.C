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

#include "pdCyclicWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdCyclicWallPatch, 0);

addToRunTimeSelectionTable(pdCyclicBoundary, pdCyclicWallPatch, dictionary);


void pdCyclicWallPatch::readProperties()
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdCyclicWallPatch::pdCyclicWallPatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdCyclicWallPatch::~pdCyclicWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void pdCyclicWallPatch::calculateProperties()
{

}

void pdCyclicWallPatch::initialConfiguration()
{}



void pdCyclicWallPatch::controlMol
(
    pdParcel& p,
    pdParcel::trackingData& td
)
{
    //- may not work in parallel, should trip switch processor before this. need to check

    const label& faceI = p.face();          //- face ID

    vector nF = mesh_.faceAreas()[faceI];   //- face vector
    nF /= mag(nF);                          //- face unit vector

    //-figure out which cyclic face we're at
    label fA = findIndex(coupledFacesA_, faceI);
    label fB = findIndex(coupledFacesB_, faceI);

    /*Info << "faces: " << faces_ << endl;
    Info << "   coupledA: " << coupledFacesA_ << endl;
    Info << "   coupledB: " << coupledFacesB_ << endl;
    Info << "   fA" << fA << endl;
    Info << "   fB" << fB << endl;*/

    //index that wasn't found is the destination patch
    vector faceD = vector::zero;
    if(fA < 0)
    {
        vector faceD = mesh_.faceCentres()[coupledFacesA_[fB]];
    }
    else
    {
        vector faceD = mesh_.faceCentres()[coupledFacesB_[fA]];
    }

    /*Info << "position: " << p.position() << endl;
    Info << "    faceI: " << faceI << endl;
    Info << "    nF: " << nF << endl;
    Info << "    faceD: " << faceD << endl;*/

    if (nF.x() != 0.0)
    {
        //Info << "faceD.x(): " << faceD.x() << endl;
        //Info << "pos.x()" << p.position().x() << endl;
        //p.position().x() = faceD.x();

        if (p.position().x() <= 0.0)
        {
            p.position().x() = 6.283;
            //p.position().x() = faceD.x();
        }
        else if (p.position().x() >= 6.283)
        {
            p.position().x() = 0.0;
        }

    }
    else if (nF.y() != 0.0)
    {
        //Info << "faceD.y(): " << faceD.y() << endl;
        //Info << "pos.y()" << p.position().y() << endl;
        //p.position().y() = faceD.x();

        if (p.position().y() <= 0.0)
        {
            p.position().y() = 0.011;
        }
        else if (p.position().y() >= 0.011)
        {
            p.position().y() = 0.0;
        }
    }
    else if (nF.z() != 0.0)
    {
        p.position().z() = faceD.z();
    }
    /*Info << "new position: " << p.position() << endl;
    Info << "old velocity: " << p.U() << endl;*/
}



void pdCyclicWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void pdCyclicWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    readProperties();
}




} // End namespace Foam

// ************************************************************************* //
