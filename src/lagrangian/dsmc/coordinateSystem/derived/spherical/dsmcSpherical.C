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
    dsmcSpherical

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcSpherical.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcSpherical, 0);

    addToRunTimeSelectionTable
    (
        dsmcCoordinateSystem,
        dsmcSpherical,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcSpherical::sphericalWeighting()
{
    forAll(cloud_.cellOccupancy(), c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cloud_.cellOccupancy()[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel& p = *molsInCell[mIC];

            const scalar oldRadialWeight = p.RWF();

            const scalar newRadialWeight = RWF(c);

            const scalar rwfFactor = newRadialWeight / oldRadialWeight;

            p.RWF() = newRadialWeight;

            if (oldRadialWeight > newRadialWeight)
            {
                //- particle might be cloned
                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;

                while(prob > 1.0)
                {
                    // Add a particle and reduce the probability by 1. If a
                    // stuck (adsorbed) parcel is cloned, the new parcel(s)
                    // should also be stuck.
                    if (p.isFree())
                    {
                        vector U = p.U();

                        cloud_.addNewParcel
                        (
                            p.position(),
                            U,
                            p.RWF(),
                            p.ERot(),
                            p.ELevel(),
                            p.cell(),
                            p.tetFace(),
                            p.tetPt(),
                            p.typeId(),
                            p.newParcel(),
                            p.classification(),
                            p.vibLevel()
                        );

                        prob -= 1.0;
                    }
                    else
                    {
                        // this particle is stuck -> cloned particles will also
                        // be initialized as stuck with the wall properties of
                        // the parent parcel:

                        scalarField wallTemperature =
                            p.stuck().wallTemperature();
                        vectorField wallVectors =
                            p.stuck().wallVectors();

                        // adjust the pre-interaction energy and momentum to
                        // the new RWF, factor is [0, 1]:
                        wallTemperature[3] *= rwfFactor;
                        wallVectors[3] *= rwfFactor;

                        cloud_.addNewStuckParcel
                        (
                            p.position(),
                            vector::zero,
                            p.RWF(),
                            p.ERot(),
                            p.ELevel(),
                            p.cell(),
                            p.tetFace(),
                            p.tetPt(),
                            p.typeId(),
                            p.newParcel(),
                            p.classification(),
                            p.vibLevel(),
                            wallTemperature,
                            wallVectors
                        );

                        prob -= 1.0;
                    }
                }
                // Add another parcel with probability equal to the remainder
                if (prob > cloud_.rndGen().sample01<scalar>())
                {
                    // If a stuck (adsorbed) parcel is cloned, the new parcel(s)
                    // should also be stuck.
                    if (p.isFree())
                    {
                        vector U = p.U();

                        cloud_.addNewParcel
                        (
                            p.position(),
                            U,
                            p.RWF(),
                            p.ERot(),
                            p.ELevel(),
                            p.cell(),
                            p.tetFace(),
                            p.tetPt(),
                            p.typeId(),
                            p.newParcel(),
                            p.classification(),
                            p.vibLevel()
                        );
                    }
                    else
                    {
                        // this particle is stuck -> cloned particles will also
                        // be initialized as stuck with the wall properties of
                        // the parent parcel:

                        scalarField wallTemperature =
                            p.stuck().wallTemperature();
                        vectorField wallVectors =
                            p.stuck().wallVectors();

                        // adjust the pre-interaction energy and momentum to
                        // the new RWF, factor is [0, 1]:
                        wallTemperature[3] *= rwfFactor;
                        wallVectors[3] *= rwfFactor;

                        cloud_.addNewStuckParcel
                        (
                            p.position(),
                            vector::zero,
                            p.RWF(),
                            p.ERot(),
                            p.ELevel(),
                            p.cell(),
                            p.tetFace(),
                            p.tetPt(),
                            p.typeId(),
                            p.newParcel(),
                            p.classification(),
                            p.vibLevel(),
                            wallTemperature,
                            wallVectors
                        );
                    }
                }
                // if the original parcel was stuck we also have to adjust its
                // wall properties:
                if (p.isStuck())
                {
                    p.stuck().wallTemperature()[3] *= rwfFactor;
                    p.stuck().wallVectors()[3] *= rwfFactor;
                }
            }
            else if (newRadialWeight > oldRadialWeight)
            {
                //- particle might be deleted
                if ((oldRadialWeight/newRadialWeight) < cloud_.rndGen().sample01<scalar>())
                {
                    cloud_.deleteParticle(p);
                }
            }
        }
    }
}


void dsmcSpherical::updateRWF()
{
    forAll(RWF_, celli)
    {
        RWF_[celli] = recalculateRWF(celli);
    }

    forAll(RWF_.boundaryField(), patchi)
    {
        fvPatchScalarField& pRWF = RWF_.boundaryFieldRef()[patchi];

        forAll(pRWF, facei)
        {
            const label celli =
                mesh_.boundaryMesh()[patchi].faceCells()[facei];

            pRWF[facei] = RWF_[celli];
        }
    }
}


scalar dsmcSpherical::recalculateRWF(const label cellI) const
{
    scalar RWF = 1.0;

    if (rWMethod_ == "particleAverage")
    {
        const DynamicList<dsmcParcel*>& cellParcels(cloud_.cellOccupancy()[cellI]);

        RWF = 0.0;
        label nMols = 0;

        forAll(cellParcels, i)
        {
            const dsmcParcel& p = *cellParcels[i];

            const scalar radius =
                sqrt
                (
                    sqr(p.position().x() - origin_.x())
                  + sqr(p.position().y() - origin_.y())
                  + sqr(p.position().z() - origin_.z())
                );

            RWF += 1.0 + (maxRWF() - 1.0)*sqr(radius/radialExtent());

            nMols++;
        }

        RWF /= max(nMols, 1);
    }
    else
    {
        const point& cC = mesh_.cellCentres()[cellI];
        const scalar radius =
            sqrt
            (
                sqr(cC.x() - origin_.x())
              + sqr(cC.y() - origin_.y())
              + sqr(cC.z() - origin_.z())
            );

        RWF += (maxRWF() - 1.0)*sqr(radius/radialExtent());
    }

    return RWF;
}


void dsmcSpherical::writeSphericalInfo() const
{
    Info<< nl << "Spherical simulation:" << nl
        << "- coordinate system origin" << tab << origin_ << nl
        << "- radial weighting method" << tab << rWMethod_ << "-based" << nl
        << "- radial extent" << tab << radialExtent_ << nl
        << "- maximum radial weighting factor" << tab << maxRWF_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcSpherical::dsmcSpherical
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcCoordinateSystem(t, mesh, cloud),
    cloud_(cloud),
    radialExtent_(0.0),
    maxRWF_(1.0),
    origin_(vector::zero),
    rWMethod_(word::null),
    RWF_
    (
        IOobject
        (
            "RWF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("RWF", dimless, 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSpherical::~dsmcSpherical()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcSpherical::checkCoordinateSystemInputs(const bool init)
{
    rWMethod_ = cloud_.particleProperties().subDict("sphericalProperties")
        .lookupOrDefault<word>("radialWeightingMethod", "cell");

    if (rWMethod_ != "cell" and rWMethod_ != "particleAverage")
    {
        FatalErrorInFunction
            << "The radial weighting method is badly defined. Choices in "
            << "constant/dsmcProperties are \"cell\" or \"particleAverage\". "
            << "Please edit the entry: radialWeightingMethod."
            << exit(FatalError);
    }

    maxRWF_ = readScalar
        (
            cloud_.particleProperties().subDict("sphericalProperties")
                .lookup("maxRadialWeightingFactor")
        );

    origin_ = cloud_.particleProperties().subDict("sphericalProperties")
        .lookupOrDefault<vector>("origin", vector::zero);

    scalarField radii(mesh_.faceCentres().size(), 0.0);

    forAll(mesh_.faceCentres(), i)
    {
        radii[i] = sqrt
            (
                sqr(mesh_.faceCentres()[i].x() - origin_.x())
              + sqr(mesh_.faceCentres()[i].y() - origin_.y())
              + sqr(mesh_.faceCentres()[i].z() - origin_.z())
            );
    }

    radialExtent_ = gMax(radii);

    writeCoordinateSystemInfo();

    // The RWFs have to be initialized only during the first run. This is done
    // by dsmcInitialise+. After that the RWFs will be read from file during
    // construction (this is for example the case in repeated running of
    // dsmcFoam+ when using dynamic load balancing).
    if (init)
    {
        // "particleAverage" cannot be used in dsmcInitialise+, "cell" is thus
        // employed
        rWMethod_ = "cell";

        updateRWF();
    }
}


void dsmcSpherical::evolve()
{
    if (rWMethod_ == "particleAverage")
    {
        updateRWF();
    }

    sphericalWeighting();
    cloud_.reBuildCellOccupancy();
}


void dsmcSpherical::writeCoordinateSystemInfo() const
{
    writeSphericalInfo();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
