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
{}


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
    rhoNIntBF_(),
    rhoNElecBF_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    mccSpeciesBF_(),
    momentumBF_(),
    UMeanBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    qBF_(),
    fDBF_(),
    evmsBF_(),
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
{}


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
        nParticlesOnStickingBoundaries_[patchi][facei] = pnStuckParticles[facei];
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
    if(mesh_.time().outputTime())
    {
        if(cloud_.boundaries().isAStickingPatch())
        {
            writenStuckParticles();
        }

        if(cloud_.boundaries().isAAbsorbingPatch())
        {
            writenAbsorbedParticles();
        }

        if(cloud_.boundaries().isAFieldPatch())
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
        forAll(rhoNBF_[i], j)
        {
            rhoNIntBF_[i][j] = 0.0;
            rhoNElecBF_[i][j] = 0.0;
            rhoNBF_[i][j] = 0.0;
            rhoMBF_[i][j] = 0.0;
            linearKEBF_[i][j] = 0.0;
            mccSpeciesBF_[i][j] = 0.0;
            momentumBF_[i][j] = vector::zero;
            UMeanBF_[i][j] = vector::zero;
            rotationalEBF_[i][j] = 0.0;
            rotationalDofBF_[i][j] = 0.0;
            vibrationalEBF_[i][j] = 0.0;
            electronicEBF_[i][j] = 0.0;
            qBF_[i][j] = 0.0;
            fDBF_[i][j] = vector::zero;
        }

        forAll(evmsBF_[i], m)
        {
            forAll(evmsBF_[i][m], j)
            {
                evmsBF_[i][m][j] = 0.0;
            }
        }
    }

    forAll(nParticlesOnStickingBoundaries_, patchi)
    {
        forAll(nParticlesOnStickingBoundaries_[patchi], facei)
        {
            nParticlesOnStickingBoundaries_[patchi][facei] = 0.0;
        }
    }
}


void boundaryMeasurements::reset()
{
    //- reset sizes of the fields after mesh is changed
    const label nSpecies = typeIds_.size();
    const label nPatches = mesh_.boundaryMesh().size();

    rhoNIntBF_.setSize(nSpecies);
    rhoNElecBF_.setSize(nSpecies);
    rhoNBF_.setSize(nSpecies);
    rhoMBF_.setSize(nSpecies);
    linearKEBF_.setSize(nSpecies);
    mccSpeciesBF_.setSize(nSpecies);
    momentumBF_.setSize(nSpecies);
    UMeanBF_.setSize(nSpecies);
    rotationalEBF_.setSize(nSpecies);
    rotationalDofBF_.setSize(nSpecies);
    vibrationalEBF_.setSize(nSpecies);
    electronicEBF_.setSize(nSpecies);
    qBF_.setSize(nSpecies);
    fDBF_.setSize(nSpecies);
    evmsBF_.setSize(nSpecies);

    forAll(typeIds_, i)
    {
        rhoNIntBF_[i].setSize(nPatches);
        rhoNElecBF_[i].setSize(nPatches);
        rhoNBF_[i].setSize(nPatches);
        rhoMBF_[i].setSize(nPatches);
        linearKEBF_[i].setSize(nPatches);
        mccSpeciesBF_[i].setSize(nPatches);
        momentumBF_[i].setSize(nPatches);
        UMeanBF_[i].setSize(nPatches);
        rotationalEBF_[i].setSize(nPatches);
        rotationalDofBF_[i].setSize(nPatches);
        vibrationalEBF_[i].setSize(nPatches);
        electronicEBF_[i].setSize(nPatches);
        qBF_[i].setSize(nPatches);
        fDBF_[i].setSize(nPatches);

        evmsBF_[i].setSize(cloud_.constProps(typeIds_[i]).thetaV().size());
        forAll(evmsBF_[i], m)
        {
            evmsBF_[i][m].setSize(nPatches);
        }

        forAll(rhoNBF_[i], j)
        {
            const label nFaces = mesh_.boundaryMesh()[j].size();

            rhoNIntBF_[i][j].setSize(nFaces, 0.0);
            rhoNElecBF_[i][j].setSize(nFaces, 0.0);
            rhoNBF_[i][j].setSize(nFaces, 0.0);
            rhoMBF_[i][j].setSize(nFaces, 0.0);
            linearKEBF_[i][j].setSize(nFaces, 0.0);
            mccSpeciesBF_[i][j].setSize(nFaces, 0.0);
            momentumBF_[i][j].setSize(nFaces, vector::zero);
            UMeanBF_[i][j].setSize(nFaces, vector::zero);
            rotationalEBF_[i][j].setSize(nFaces, 0.0);
            rotationalDofBF_[i][j].setSize(nFaces, 0.0);
            vibrationalEBF_[i][j].setSize(nFaces, 0.0);
            electronicEBF_[i][j].setSize(nFaces, 0.0);
            qBF_[i][j].setSize(nFaces, 0.0);
            fDBF_[i][j].setSize(nFaces, vector::zero);
        }

        forAll(evmsBF_[i], m)
        {
            forAll(evmsBF_[i][m], j)
            {
                const label nFaces = mesh_.boundaryMesh()[j].size();
                evmsBF_[i][m][j].setSize(nFaces, 0.0);
            }
        }
    }

    forAll(nParticlesOnStickingBoundaries_, patchi)
    {
        forAll(nParticlesOnStickingBoundaries_[patchi], facei)
        {
            nParticlesOnStickingBoundaries_[patchi][facei] = 0.0;
            nAbsorbedParticles_[patchi][facei] = 0.0;
            boundaryT_[patchi][facei] = 0.0;
            boundaryU_[patchi][facei] = vector::zero;
        }
    }
}


void boundaryMeasurements::updateFields(dsmcParcel& p)
{}


}  // End namespace Foam

// ************************************************************************* //
