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
#include "dsmcCloud.H"

namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void boundaryMeasurements::writenStuckParcels()
{
    tmp<volScalarField> tnStuckParcels
    (
        new volScalarField
        (
            IOobject
            (
                "nStuckParcels",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nStuckParcels", dimless, 0.0)
        )
    );
    
    volScalarField& nStuckParcels = tnStuckParcels.ref();
    
    nStuckParcels.boundaryFieldRef() = nParcelsOnStickingBoundaries_;
    
    nStuckParcels.write();
}


void boundaryMeasurements::writenAbsorbedParcels()
{
    tmp<volScalarField> tnAbsorbedParcels
    (
        new volScalarField
        (
            IOobject
            (
                "nAbsorbedParcels",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nAbsorbedParcels", dimless, 0.0)
        )
    );
    
    volScalarField& nAbsorbedParcels = tnAbsorbedParcels.ref();
    
    nAbsorbedParcels.boundaryFieldRef() = nAbsorbedParcels_;
    
    nAbsorbedParcels.write();
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
    nParcelsOnStickingBoundaries_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    nAbsorbedParcels_
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
    nParcelsOnStickingBoundaries_
    (
        mesh_.boundary(),
        mesh_.V(),
        calculatedFvPatchScalarField::typeName
    ),
    nAbsorbedParcels_
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

void boundaryMeasurements::updatenStuckParcelOnPatch
(
    const label patchi,
    const scalarList& pnStuckParcels
)
{
    forAll(pnStuckParcels, facei)
    {
        nParcelsOnStickingBoundaries_[patchi][facei] = 
            pnStuckParcels[facei];
    }    
}


void boundaryMeasurements::updatenAbsorbedParcelOnPatch
(
    const label patchi,
    const label facei
)
{
    nAbsorbedParcels_[patchi][facei]++;
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


void boundaryMeasurements::setBoundarynStuckParcels
(
    const label patchi,
    const scalarList& pnStuckParcels
)
{
    forAll(pnStuckParcels, facei)
    {
        nParcelsOnStickingBoundaries_[patchi][facei] = pnStuckParcels[facei];
    }    
}


void boundaryMeasurements::setBoundarynAbsorbedParcels
(
    const label patchi,
    const scalarList& pnAbsorbedParcels
)
{
    forAll(pnAbsorbedParcels, facei)
    {
        nAbsorbedParcels_[patchi][facei] = pnAbsorbedParcels[facei];
    }    
}


void boundaryMeasurements::outputResults()
{
    if(mesh_.time().outputTime())
    {
        if(cloud_.boundaries().isAStickingPatch())
        {
            writenStuckParcels();
        }
        
        if(cloud_.boundaries().isAAbsorbingPatch())
        {
            writenAbsorbedParcels();
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
    
    forAll(nParcelsOnStickingBoundaries_, patchi)
    {
        forAll(nParcelsOnStickingBoundaries_[patchi], facei)
        {
            nParcelsOnStickingBoundaries_[patchi][facei] = 0.0;
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
    
    forAll(nParcelsOnStickingBoundaries_, patchi)
    {
        forAll(nParcelsOnStickingBoundaries_[patchi], facei)
        {
            nParcelsOnStickingBoundaries_[patchi][facei] = 0.0;
            nAbsorbedParcels_[patchi][facei] = 0.0;
            boundaryT_[patchi][facei] = 0.0;
            boundaryU_[patchi][facei] = vector::zero;
        }
    }
}


void boundaryMeasurements::updateFields(dsmcParcel& p)
{}


}  // End namespace Foam

// ************************************************************************* //
