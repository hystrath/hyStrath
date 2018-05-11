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
            << "Cannot have zero typeIds being adsorbed." << nl << "in: "
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
            saturationLimitPerSquareMeters*facesArea[facei]
           /cloud_.nParticles(patchId(), facei);
    }
}


void dsmcAbsorbingWallPatch::readPatchField()
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
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nAbsorbedParcels", dimless, 0.0)
        )
    );
    
    volScalarField& nAbsorbedParcels = tnAbsorbedParcels.ref();
    
    cloud_.boundaryFluxMeasurements().setBoundarynAbsorbedParcels
    (
        patchId(),
        nAbsorbedParcels.boundaryField()[patchId()]
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
            .pnAbsorbedParcels(patchi)[facei] > 1e-6
    )
    {
        return true;
    }
    
    return false;
}


void dsmcAbsorbingWallPatch::absorbParticle
(
    const label patchi,
    const label facei,
    dsmcParcel::trackingData& td
)
{
    //- Delete particle
    td.keepParticle = false;
    
    //- Update the boundaryMeasurement relative to this absorbing patch
    cloud_.boundaryFluxMeasurements().updatenAbsorbedParcelOnPatch
    (
        patchi,
        facei
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
        const scalar absorptionProbability = absorptionProb(iD);
            
        if
        (
            absorptionProbability > cloud_.rndGen().scalar01()
         && isNotSaturated(wppIndex, wppLocalFace)
        )
        {
            //- absorb particle
            absorbParticle(wppIndex, wppLocalFace, td);
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
