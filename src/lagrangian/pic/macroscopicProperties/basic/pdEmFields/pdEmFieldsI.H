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


Description

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Access
/*
//template <class pdParcel>
inline volScalarField::Boundary& pdEmFields::lambdaDBF()
{
    return lambdaD_.boundaryFieldRef();
}*/

//template <class pdParcel>
inline volScalarField::Boundary& pdEmFields::rhoQBF()
{
    return rhoQ_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdEmFields::wallQBF()
{
    return wallQ_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdEmFields::rhoEFBF()
{
    return rhoEF_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdEmFields::phiEBF()
{
    return phiE_.boundaryFieldRef();
}

//template <class pdParcel>
inline volVectorField::Boundary& pdEmFields::EBF()
{
    return E_.boundaryFieldRef();
}

//template <class pdParcel>
inline volVectorField::Boundary& pdEmFields::BBF()
{
    return B_.boundaryFieldRef();
}

//template <class pdParcel>
inline volVectorField::Boundary& pdEmFields::JpBF()
{
    return Jp_.boundaryFieldRef();
}
/*
//template <class pdParcel>
inline const volScalarField& pdEmFields::lambdaD() const
{
    return lambdaD_;
}*/

//template <class pdParcel>
inline const volScalarField& pdEmFields::rhoQ() const
{
    return rhoQ_;
}

//template <class pdParcel>
inline const volScalarField& pdEmFields::wallQ() const
{
    return wallQ_;
}

//template <class pdParcel>
inline const volScalarField& pdEmFields::rhoEF() const
{
    return rhoEF_;
}

//template <class pdParcel>
inline const volScalarField& pdEmFields::phiE() const
{
    return phiE_;
}

//template <class pdParcel>
inline const volVectorField& pdEmFields::E() const
{
    return E_;
}

//template <class pdParcel>
inline const volVectorField& pdEmFields::B() const
{
    return B_;
}

//template <class pdParcel>
inline const volVectorField& pdEmFields::Jp() const
{
    return Jp_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
