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

Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Access

inline scalar boundaryMeasurements::speciesRhoNIntBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesRhoNIntBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesRhoNElecBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesRhoNElecBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesRhoNBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesRhoNBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesRhoMBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesRhoMBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesLinearKEBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesLinearKEBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesMccBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesMccBF_[speciei][patchi][facei];
}


inline vector boundaryMeasurements::speciesMomentumBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesMomentumBF_[speciei][patchi][facei];
}


inline vector boundaryMeasurements::speciesUMeanBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesUMeanBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesErotBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesErotBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesZetaRotBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesZetaRotBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesEvibBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesEvibBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesEvibModBF
(
    const label speciei,
    const label modi,
    const label patchi,
    const label facei
) const
{
    return speciesEvibModBF_[speciei][modi][patchi][facei];
}


inline scalar boundaryMeasurements::speciesEelecBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesEelecBF_[speciei][patchi][facei];
}


inline scalar boundaryMeasurements::speciesqBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesqBF_[speciei][patchi][facei];
}


inline vector boundaryMeasurements::speciesfDBF
(
    const label speciei,
    const label patchi,
    const label facei
) const
{
    return speciesfDBF_[speciei][patchi][facei];
}


inline const fvPatchScalarField& boundaryMeasurements::pnStuckParticles
(
    const label patchi
) const
{
    return nParticlesOnStickingBoundaries_[patchi];
}


inline const fvPatchScalarField& boundaryMeasurements::pnAbsorbedParticles
(
    const label patchi
) const
{
    return nAbsorbedParticles_[patchi];
}


inline const fvPatchScalarField& boundaryMeasurements::boundaryT
(
    const label patchi
) const
{
    return boundaryT_[patchi];
}


inline const fvPatchVectorField& boundaryMeasurements::boundaryU
(
    const label patchi
) const
{
    return boundaryU_[patchi];
}


//- Edit

inline scalar& boundaryMeasurements::speciesRhoNIntBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesRhoNIntBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesRhoNElecBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesRhoNElecBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesRhoNBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesRhoNBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesRhoMBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesRhoMBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesLinearKEBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesLinearKEBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesMccBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesMccBF_[speciei][patchi][facei];
}


inline vector& boundaryMeasurements::speciesMomentumBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesMomentumBF_[speciei][patchi][facei];
}


inline vector& boundaryMeasurements::speciesUMeanBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesUMeanBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesErotBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesErotBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesZetaRotBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesZetaRotBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesEvibBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesEvibBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesEvibModBF
(
    const label speciei,
    const label modi,
    const label patchi,
    const label facei
)
{
    return speciesEvibModBF_[speciei][modi][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesEelecBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesEelecBF_[speciei][patchi][facei];
}


inline scalar& boundaryMeasurements::speciesqBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesqBF_[speciei][patchi][facei];
}


inline vector& boundaryMeasurements::speciesfDBF
(
    const label speciei,
    const label patchi,
    const label facei
)
{
    return speciesfDBF_[speciei][patchi][facei];
}


inline fvPatchScalarField& boundaryMeasurements::pnStuckParticles
(
    const label patchi
)
{
    return nParticlesOnStickingBoundaries_[patchi];
}


inline fvPatchScalarField& boundaryMeasurements::pnAbsorbedParticles
(
    const label patchi
)
{
    return nAbsorbedParticles_[patchi];
}


inline fvPatchScalarField& boundaryMeasurements::boundaryT
(
    const label patchi
)
{
    return boundaryT_[patchi];
}


inline fvPatchVectorField& boundaryMeasurements::boundaryU
(
    const label patchi
)
{
    return boundaryU_[patchi];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
