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

#include "dsmcAbsorbingWallPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAbsorbingWallPatch, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcAbsorbingWallPatch::setProperties()
{
    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcAbsorbingWallPatch::setProperties()")
            << "Cannot have zero typeIds being absorbed." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        const label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcAbsorbingWallPatch::setProperties()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    absorptionProbs_.clear();

    absorptionProbs_.setSize(typeIds_.size(), 0.0);

    forAll(absorptionProbs_, i)
    {
        absorptionProbs_[i] = readScalar
        (
            propsDict_.subDict("absorptionProbabilities")
                .lookup(moleculesReduced[i])
        );
    }

    const scalar saturationLimitPerSquareMeters = propsDict_
        .lookupOrDefault<scalar>("saturationLimit", VGREAT);

    const scalarList& facesArea = mesh_.magSf().boundaryField()[patchId()];

    forAll(saturationLimit_, facei)
    {
        saturationLimit_[facei] =
            saturationLimitPerSquareMeters*facesArea[facei];
    }
}


void dsmcAbsorbingWallPatch::readPatchField()
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
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nAbsorbedParticles", dimless, 0.0)
        )
    );

    volScalarField& nAbsorbedParticles = tnAbsorbedParticles.ref();

    cloud_.boundaryFluxMeasurements().setBoundarynAbsorbedParticles
    (
        patchId(),
        nAbsorbedParticles.boundaryField()[patchId()]
    );
}


bool dsmcAbsorbingWallPatch::isNotSaturated
(
    const label patchi,
    const label facei
)
{
    if
    (
        saturationLimit(facei) - cloud_.boundaryFluxMeasurements()
            .pnAbsorbedParticles(patchi)[facei] > 1e-6
    )
    {
        return true;
    }

    return false;
}


void dsmcAbsorbingWallPatch::absorbParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    // Particle location on patch / face
    const label wppIndex = patchId();
    const label wppLocalFace =
        mesh_.boundaryMesh()[wppIndex].whichFace(p.face());

    // Calculate the real number of particles that this parcel represented
    const scalar nAbsorbedParticles = p.RWF() * cloud_.coordSystem().dtModel()
        .nParticles(wppIndex, wppLocalFace);

    // Delete parcel
    td.keepParticle = false;

    // Update the boundaryMeasurement relative to this absorbing patch
    cloud_.boundaryFluxMeasurements().updatenAbsorbedParticlesOnPatch
    (
        wppIndex,
        wppLocalFace,
        nAbsorbedParticles
    );
}


void dsmcAbsorbingWallPatch::testParticleForAbsorption
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    const label iD = findIndex(dsmcAbsorbingWallPatch::typeIds_, p.typeId());

    if(iD != -1)
    {
        const label wppIndex = patchId();

        const label wppLocalFace =
            mesh_.boundaryMesh()[wppIndex].whichFace(p.face());

        //- absorption probability
        const scalar absorptionProbability = absorptionProbs_[iD];

        if
        (
            absorptionProbability > cloud_.rndGen().sample01<scalar>()
         && isNotSaturated(wppIndex, wppLocalFace)
        )
        {
            //- absorb particle
            absorbParticle(p, td);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAbsorbingWallPatch::dsmcAbsorbingWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    absorptionProbs_(),
    saturationLimit_(mesh_.boundaryMesh()[patchId()].size(), 0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    readPatchField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAbsorbingWallPatch::~dsmcAbsorbingWallPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAbsorbingWallPatch::initialConfiguration()
{}


void dsmcAbsorbingWallPatch::calculateProperties()
{}


void dsmcAbsorbingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcAbsorbingWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
