/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::hTC2Model::mesh() const
{
    return mesh_;
}


inline const Foam::surfaceScalarField& Foam::hTC2Model::phi() const
{
    if (turbulencePtr_)
    {
        return turbulencePtr_->phi();
    }
    else
    {
        FatalErrorIn
        (
            "const Foam::compressible::turbulenceModel2& "
            "Foam::hTC2Model::turbulence() const "
        )   << "turbulencePtr_ is empty. Please use "
            << "hTC2Model::setTurbulence "
            << "(compressible::turbulenceModel2& )"
            << abort(FatalError);

        return turbulencePtr_->phi();
    }
}


inline const Foam::compressible::turbulenceModel2&
Foam::hTC2Model::turbulence() const
{
    if (turbulencePtr_)
    {
        return *turbulencePtr_;
    }
    else
    {
        FatalErrorIn
        (
            "const Foam::compressible::turbulenceModel2& "
            "Foam::hTC2Model::turbulence() const "
        )   << "turbulencePtr_ is empty. Please use "
            << "hTC2Model::setTurbulence "
            << "(compressible::turbulenceModel2& )"
            << abort(FatalError);

        return *turbulencePtr_;
    }
}


inline const Foam::Switch& Foam::hTC2Model::active() const
{
    return active_;
}


inline void Foam::hTC2Model::setTurbulence
(
    compressible::turbulenceModel2& turbModel
)
{
    turbulencePtr_ = &turbModel;
}


inline const Foam::dictionary& Foam::hTC2Model::coeffs() const
{
    return coeffs_;
}

// ************************************************************************* //
