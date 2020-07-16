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

#include "dsmcAbsorbingStickingDiffuseWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAbsorbingStickingDiffuseWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcAbsorbingStickingDiffuseWallFieldPatch,
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAbsorbingStickingDiffuseWallFieldPatch::
    dsmcAbsorbingStickingDiffuseWallFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseWallPatch(t, mesh, cloud, dict),
    dsmcAbsorbingWallPatch(t, mesh, cloud, dict),
    dsmcStickingWallPatch(t, mesh, cloud, dict)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    dsmcAbsorbingWallPatch::setProperties();
    dsmcStickingWallPatch::setProperties();
    dsmcDiffuseWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAbsorbingStickingDiffuseWallFieldPatch::
    ~dsmcAbsorbingStickingDiffuseWallFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAbsorbingStickingDiffuseWallFieldPatch::initialConfiguration()
{}


void dsmcAbsorbingStickingDiffuseWallFieldPatch::calculateProperties()
{}


void dsmcAbsorbingStickingDiffuseWallFieldPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    if(p.isFree())
    {
        measurePropertiesBeforeControl(p);

        const label iDab =
            findIndex(dsmcAbsorbingWallPatch::typeIds_, p.typeId());

        const label iDst =
            findIndex(dsmcStickingWallPatch::typeIds_, p.typeId());

        const label wppIndex = patchId();

        const label wppLocalFace =
            mesh_.boundaryMesh()[wppIndex].whichFace(p.face());

        //- Calculation of the local patch temperature
        const scalar localPatchTemperature =
            dsmcFieldPatchBoundary::patchLocalTemperature(p);

        //- Calculation of the local patch velocity
        const vector& localPatchVelocity =
            dsmcFieldPatchBoundary::patchLocalVelocity(p);

        bool performDiffusiveReflection = false;

        if(iDab != -1 && iDst != -1)
        {
            //- absorption probability
            const scalar absorptionProbability = absorptionProbs_[iDab];

            //- adsorption probability
            const scalar adsorptionProbability = adsorptionProbs_[iDst];

            if
            (
                dsmcStickingWallPatch::isNotSaturated(wppLocalFace)
             && dsmcAbsorbingWallPatch::isNotSaturated(wppIndex, wppLocalFace)
            )
            {
                const scalar sumAbsAdsProbabilities = absorptionProbability
                    + adsorptionProbability;

                if(sumAbsAdsProbabilities > cloud_.rndGen().sample01<scalar>())
                {
                    //- Either absorption or adsorption must be operated.
                    //  The probability of adsorption is rescaled
                    if
                    (
                        adsorptionProbability/sumAbsAdsProbabilities
                            > cloud_.rndGen().sample01<scalar>()
                    )
                    {
                        //- particle is adsorbed
                        dsmcStickingWallPatch::adsorbParticle
                        (
                            p,
                            getLocalTemperature(p.position()[depthAxis()])
                        );
                    }
                    else
                    {
                        //- absorb particle
                        dsmcAbsorbingWallPatch::absorbParticle(p, td);
                    }
                }
                else
                {
                    performDiffusiveReflection = true;
                }
            }
            else
            {
                //- Either one or both bounding mechanism are saturated
                if(dsmcStickingWallPatch::isNotSaturated(wppLocalFace))
                {
                    if(adsorptionProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //- particle is adsorbed
                        dsmcStickingWallPatch::adsorbParticle
                        (
                            p,
                            getLocalTemperature(p.position()[depthAxis()])
                        );
                    }
                    else
                    {
                        performDiffusiveReflection = true;
                    }
                }
                else if
                (
                    dsmcAbsorbingWallPatch::isNotSaturated
                    (
                        wppIndex,
                        wppLocalFace
                    )
                )
                {
                    if(absorptionProbability > cloud_.rndGen().sample01<scalar>())
                    {
                        //- absorb particle
                        dsmcAbsorbingWallPatch::absorbParticle(p, td);
                    }
                    else
                    {
                        performDiffusiveReflection = true;
                    }
                }
                else
                {
                    performDiffusiveReflection = true;
                }
            }
        }
        else
        {
            performDiffusiveReflection = true;
        }

        if(performDiffusiveReflection)
        {
            //- diffuse reflection
            dsmcDiffuseWallPatch::performDiffuseReflection
            (
                p,
                localPatchTemperature,
                localPatchVelocity
            );

            measurePropertiesAfterControl(p);
        }
    }

    //- Separate loop as the particle may have been stuck in the previous loop
    //  If the particle is stuck, consider parcel for release, i.e., desorption
    if(p.isStuck())
    {
        dsmcStickingWallPatch::testForDesorption(p);
    }

    //- Update the boundaryMeasurement relative to this sticking patch
    cloud_.boundaryFluxMeasurements().updatenStuckParticlesOnPatch
    (
        patchId(),
        nStuckParticles_
    );
}


void dsmcAbsorbingStickingDiffuseWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcAbsorbingStickingDiffuseWallFieldPatch::updateProperties
(
    const dictionary& newDict
)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    dsmcAbsorbingWallPatch::setProperties();
    dsmcStickingWallPatch::setProperties();
    dsmcDiffuseWallPatch::setProperties();
}


} // End namespace Foam

// ************************************************************************* //
