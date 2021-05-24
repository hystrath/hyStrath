/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

Class
    boundaryMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "boundaryMeasurements.H"

#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"

namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void boundaryMeasurements::writenStuckParticles()
{
    tmp<volScalarField> tnStuckParticles
    (
        new volScalarField
        (
            IOobject
            (
                "nStuckParticles",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nStuckParticles", dimless, 0.0)
        )
    );

    volScalarField& nStuckParticles = tnStuckParticles.ref();

    nStuckParticles.boundaryFieldRef() = nParticlesOnStickingBoundaries_;

    nStuckParticles.write();
}


void boundaryMeasurements::writenAbsorbedParticles()
{
    tmp<volScalarField> tnAbsorbedParticles
    (
        new volScalarField
        (
            IOobject
            (
                "nAbsorbedParticles",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nAbsorbedParticles", dimless, 0.0)
        )
    );

    volScalarField& nAbsorbedParticles = tnAbsorbedParticles.ref();

    nAbsorbedParticles.boundaryFieldRef() = nAbsorbedParticles_;

    nAbsorbedParticles.write();
}


void boundaryMeasurements::writePatchFields()
{
    //- Temperature field
    tmp<volScalarField> tboundaryT
    (
        new volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("boundaryT", dimTemperature, 0.0)
        )
    );

    volScalarField& boundaryT = tboundaryT.ref();

    boundaryT.boundaryFieldRef() = boundaryT_;

    boundaryT.write();

    //- Velocity field
    tmp<volVectorField> tboundaryU
    (
        new volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("boundaryU", dimVelocity, vector::zero)
        )
    );

    volVectorField& boundaryU = tboundaryU.ref();

    boundaryU.boundaryFieldRef() = boundaryU_;

    boundaryU.write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and cloud
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    nParticlesOnStickingBoundaries_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    nAbsorbedParticles_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    boundaryT_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    boundaryU_
    (
        mesh_.boundary(),
        mesh_.C(),
        calculatedFvPatchVectorField::typeName
    )
{
    nParticlesOnStickingBoundaries_ = 0.0;
    nAbsorbedParticles_ = 0.0;
}


//- Construct from mesh, cloud and boolean (dsmcFoam)
boundaryMeasurements::boundaryMeasurements
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const bool& dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(cloud_.typeIdList().size(), -1),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    speciesRhoNBF_(),
    speciesRhoMBF_(),
    speciesLinearKEBF_(),
    speciesMccBF_(),
    speciesMomentumBF_(),
    speciesUMeanBF_(),
    speciesErotBF_(),
    speciesZetaRotBF_(),
    speciesEvibBF_(),
    speciesEelecBF_(),
    speciesqBF_(),
    speciesfDBF_(),
    speciesEvibModBF_(),
    nParticlesOnStickingBoundaries_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    nAbsorbedParticles_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    boundaryT_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    boundaryU_
    (
        mesh_.boundary(),
        mesh_.C(),
        calculatedFvPatchVectorField::typeName
    )
{
    nParticlesOnStickingBoundaries_ = 0.0;
    nAbsorbedParticles_ = 0.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

boundaryMeasurements::~boundaryMeasurements()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void boundaryMeasurements::updatenStuckParticlesOnPatch
(
    const label patchi,
    const scalarList& pnStuckParticles
)
{
    forAll(pnStuckParticles, facei)
    {
        nParticlesOnStickingBoundaries_[patchi][facei] =
            pnStuckParticles[facei];
    }
}


void boundaryMeasurements::updatenAbsorbedParticlesOnPatch
(
    const label patchi,
    const label facei,
    const scalar nAbsorbedParticles
)
{
    nAbsorbedParticles_[patchi][facei] += nAbsorbedParticles;
}


void boundaryMeasurements::setBoundaryT
(
    const label patchi,
    const scalarList& pboundaryT
)
{
    forAll(pboundaryT, facei)
    {
        boundaryT_[patchi][facei] = pboundaryT[facei];
    }
}


void boundaryMeasurements::setBoundaryU
(
    const label patchi,
    const vectorList& pboundaryU
)
{
    forAll(pboundaryU, facei)
    {
        boundaryU_[patchi][facei] = pboundaryU[facei];
    }
}


void boundaryMeasurements::setBoundarynStuckParticles
(
    const label patchi,
    const scalarList& pnStuckParticles
)
{
    forAll(pnStuckParticles, facei)
    {
        nParticlesOnStickingBoundaries_[patchi][facei] =
            pnStuckParticles[facei];
    }
}


void boundaryMeasurements::setBoundarynAbsorbedParticles
(
    const label patchi,
    const scalarList& pnAbsorbedParticles
)
{
    forAll(pnAbsorbedParticles, facei)
    {
        nAbsorbedParticles_[patchi][facei] = pnAbsorbedParticles[facei];
    }
}


void boundaryMeasurements::outputResults()
{
    if (mesh_.time().outputTime())
    {
        if (cloud_.boundaries().isAStickingPatch())
        {
            writenStuckParticles();
        }

        if (cloud_.boundaries().isAAbsorbingPatch())
        {
            writenAbsorbedParticles();
        }

        if (cloud_.boundaries().isAFieldPatch())
        {
            writePatchFields();
        }
    }
}


void boundaryMeasurements::setInitialConfig()
{
    forAll(typeIds_, i)
    {
        typeIds_[i] = i;
    }

    reset();
}


void boundaryMeasurements::clean()
{
    //- clean geometric fields
    forAll(typeIds_, i)
    {
        forAll(speciesRhoNBF_[i], j)
        {
            speciesRhoNIntBF_[i][j] = 0.0;
            speciesRhoNElecBF_[i][j] = 0.0;
            speciesRhoNBF_[i][j] = 0.0;
            speciesRhoMBF_[i][j] = 0.0;
            speciesLinearKEBF_[i][j] = 0.0;
            speciesMccBF_[i][j] = 0.0;
            speciesMomentumBF_[i][j] = vector::zero;
            speciesUMeanBF_[i][j] = vector::zero;
            speciesErotBF_[i][j] = 0.0;
            speciesZetaRotBF_[i][j] = 0.0;
            speciesEvibBF_[i][j] = 0.0;
            speciesEelecBF_[i][j] = 0.0;
            speciesqBF_[i][j] = 0.0;
            speciesfDBF_[i][j] = vector::zero;
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            forAll(speciesEvibModBF_[i][mod], j)
            {
                speciesEvibModBF_[i][mod][j] = 0.0;
            }
        }
    }

    forAll(nParticlesOnStickingBoundaries_, j)
    {
        nParticlesOnStickingBoundaries_[j] = 0.0;
    }
}


void boundaryMeasurements::reset()
{
    //- reset sizes of the fields after mesh is changed
    const label nSpecies = typeIds_.size();
    const label nPatches = mesh_.boundaryMesh().size();

    speciesRhoNIntBF_.setSize(nSpecies);
    speciesRhoNElecBF_.setSize(nSpecies);
    speciesRhoNBF_.setSize(nSpecies);
    speciesRhoMBF_.setSize(nSpecies);
    speciesLinearKEBF_.setSize(nSpecies);
    speciesMccBF_.setSize(nSpecies);
    speciesMomentumBF_.setSize(nSpecies);
    speciesUMeanBF_.setSize(nSpecies);
    speciesErotBF_.setSize(nSpecies);
    speciesZetaRotBF_.setSize(nSpecies);
    speciesEvibBF_.setSize(nSpecies);
    speciesEelecBF_.setSize(nSpecies);
    speciesqBF_.setSize(nSpecies);
    speciesfDBF_.setSize(nSpecies);
    speciesEvibModBF_.setSize(nSpecies);
    
    forAll(typeIds_, i)
    {
        const label spId = typeIds_[i];
        const label nVibMod = cloud_.constProps(spId).nVibrationalModes();
        
        speciesRhoNIntBF_[i].setSize(nPatches);
        speciesRhoNElecBF_[i].setSize(nPatches);
        speciesRhoNBF_[i].setSize(nPatches);
        speciesRhoMBF_[i].setSize(nPatches);
        speciesLinearKEBF_[i].setSize(nPatches);
        speciesMccBF_[i].setSize(nPatches);
        speciesMomentumBF_[i].setSize(nPatches);
        speciesUMeanBF_[i].setSize(nPatches);
        speciesErotBF_[i].setSize(nPatches);
        speciesZetaRotBF_[i].setSize(nPatches);
        speciesEvibBF_[i].setSize(nPatches);
        speciesEelecBF_[i].setSize(nPatches);
        speciesqBF_[i].setSize(nPatches);
        speciesfDBF_[i].setSize(nPatches);

        speciesEvibModBF_[i].setSize(nVibMod);
        forAll(speciesEvibModBF_[i], mod)
        {
            speciesEvibModBF_[i][mod].setSize(nPatches);
        }

        forAll(speciesRhoNBF_[i], j)
        {
            const label nFaces = mesh_.boundaryMesh()[j].size();

            speciesRhoNIntBF_[i][j].setSize(nFaces, 0.0);
            speciesRhoNElecBF_[i][j].setSize(nFaces, 0.0);
            speciesRhoNBF_[i][j].setSize(nFaces, 0.0);
            speciesRhoMBF_[i][j].setSize(nFaces, 0.0);
            speciesLinearKEBF_[i][j].setSize(nFaces, 0.0);
            speciesMccBF_[i][j].setSize(nFaces, 0.0);
            speciesMomentumBF_[i][j].setSize(nFaces, vector::zero);
            speciesUMeanBF_[i][j].setSize(nFaces, vector::zero);
            speciesErotBF_[i][j].setSize(nFaces, 0.0);
            speciesZetaRotBF_[i][j].setSize(nFaces, 0.0);
            speciesEvibBF_[i][j].setSize(nFaces, 0.0);
            speciesEelecBF_[i][j].setSize(nFaces, 0.0);
            speciesqBF_[i][j].setSize(nFaces, 0.0);
            speciesfDBF_[i][j].setSize(nFaces, vector::zero);
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            forAll(speciesEvibModBF_[i][mod], j)
            {
                const label nFaces = mesh_.boundaryMesh()[j].size();
                speciesEvibModBF_[i][mod][j].setSize(nFaces, 0.0);
            }
        }
    }
    
    nParticlesOnStickingBoundaries_.setSize(nPatches);
    nAbsorbedParticles_.setSize(nPatches);
    boundaryT_.setSize(nPatches);
    boundaryU_.setSize(nPatches);

    forAll(mesh_.boundaryMesh(), j)
    {
        const label nFaces = mesh_.boundaryMesh()[j].size();
        
        nParticlesOnStickingBoundaries_[j].setSize(nFaces, 0.0);
        nAbsorbedParticles_[j].setSize(nFaces, 0.0);
        boundaryT_[j].setSize(nFaces, 0.0);
        boundaryU_[j].setSize(nFaces, vector::zero);
    }
}


void boundaryMeasurements::updateFields(dsmcParcel& p)
{}


}  // End namespace Foam

// ************************************************************************* //
