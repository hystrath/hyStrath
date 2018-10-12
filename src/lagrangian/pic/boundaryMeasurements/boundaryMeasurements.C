/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    boundaryMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "boundaryMeasurements.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "pdCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from mesh and cloud
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    pdCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{
}


//- Construct from mesh, cloud and boolean (pdFoam)
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    pdCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(cloud_.typeIdList().size(),-1),
    rhoNIntBF_(),
    rhoNBF_(),
    rhoMBF_(),
    JpBF_(),
    rhoQBF_(),
    wallQBF_(),
    linearKEBF_(),
    momentumBF_(),
    UMeanBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    vibrationalEBF_(),
    qBF_(),
    fDBF_()
{

    forAll(typeIds_, i)
    {
        typeIds_[i] = i;
    }

    rhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    rhoNBF_.setSize(typeIds_.size());
    rhoMBF_.setSize(typeIds_.size());
    JpBF_.setSize(typeIds_.size());
    rhoQBF_.setSize(typeIds_.size());
    wallQBF_.setSize(typeIds_.size());
    linearKEBF_.setSize(typeIds_.size());
    momentumBF_.setSize(typeIds_.size());
    UMeanBF_.setSize(typeIds_.size());
    rotationalEBF_.setSize(typeIds_.size());
    rotationalDofBF_.setSize(typeIds_.size());
    vibrationalEBF_.setSize(typeIds_.size());
    qBF_.setSize(typeIds_.size());
    fDBF_.setSize(typeIds_.size());

    forAll(rhoNIntBF_, i)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[i];
        rhoNIntBF_[i].setSize(patch.size(),0.0);
    }

    forAll(rhoNBF_, i)
    {
        rhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoMBF_[i].setSize(mesh_.boundaryMesh().size());
        JpBF_[i].setSize(mesh_.boundaryMesh().size());
        rhoQBF_[i].setSize(mesh_.boundaryMesh().size());
        wallQBF_[i].setSize(mesh_.boundaryMesh().size());
        linearKEBF_[i].setSize(mesh_.boundaryMesh().size());
        momentumBF_[i].setSize(mesh_.boundaryMesh().size());
        UMeanBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        rotationalDofBF_[i].setSize(mesh_.boundaryMesh().size());
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        qBF_[i].setSize(mesh_.boundaryMesh().size());
        fDBF_[i].setSize(mesh_.boundaryMesh().size());
        //****//

        forAll(rhoNBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            rhoNBF_[i][j].setSize(patch.size(),0.0);
            rhoMBF_[i][j].setSize(patch.size(),0.0);
            JpBF_[i][j].setSize(patch.size(),vector::zero);
            rhoQBF_[i][j].setSize(patch.size(),0.0);
            wallQBF_[i][j].setSize(patch.size(),0.0);
            linearKEBF_[i][j].setSize(patch.size(),0.0);
            momentumBF_[i][j].setSize(patch.size(),vector::zero);
            UMeanBF_[i][j].setSize(patch.size(),vector::zero);
            rotationalEBF_[i][j].setSize(patch.size(),0.0);
            rotationalDofBF_[i][j].setSize(patch.size(),0.0);
            vibrationalEBF_[i][j].setSize(patch.size(),0.0);
            qBF_[i][j].setSize(patch.size(),0.0);
            fDBF_[i][j].setSize(patch.size(),vector::zero);
            //****//
        }
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryMeasurements::~boundaryMeasurements()
{
    //Info << "boundaryMeasurements Destructor" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void boundaryMeasurements::clean()
{
    //- clean geometric fields

    forAll(rhoNBF_, i)
    {
        forAll(rhoNBF_[i], j)
        {
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            JpBF_[i][j] = vector::zero;
            rhoQBF_[i][j] = 0.0;
            wallQBF_[i][j] = 0.0;
            linearKEBF_[i][j] = 0.0;
            momentumBF_[i][j] = vector::zero;
            UMeanBF_[i][j] = vector::zero;
            rotationalEBF_[i][j] = 0.0;
            rotationalDofBF_[i][j] = 0.0;
            vibrationalEBF_[i][j] = 0.0;
            qBF_[i][j] = 0.0;
            fDBF_[i][j] = vector::zero;
            //****//
            rhoNIntBF_[j] = 0.0;
        }
    }
}


void boundaryMeasurements::updateFields
(
    pdParcel& p
)
{
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
