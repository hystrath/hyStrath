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


//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::qBF()
{
    return q_.boundaryFieldRef();
}


//template <class pdParcel>
inline volVectorField::Boundary& pdStandardFields::fDBF()
{
    return fD_.boundaryFieldRef();
}


//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::rhoNBF()
{
    return rhoN_.boundaryFieldRef();
}


//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::rhoMBF()
{
    return rhoM_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::linearKEBF()
{
    return linearKE_.boundaryFieldRef();
}


//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::rotationalEBF()
{
    return rotationalE_.boundaryFieldRef();
}


//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::rotationalDofBF()
{
    return rotationalDof_.boundaryFieldRef();
}

inline volScalarField::Boundary& pdStandardFields::vibrationalEBF()
{
    return vibrationalE_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::vibrationalDofBF()
{
    return vibrationalDof_.boundaryFieldRef();
}

//template <class pdParcel>
inline volVectorField::Boundary& pdStandardFields::momentumBF()
{
    return momentum_.boundaryFieldRef();
}
/*
//template <class pdParcel>
inline volVectorField::Boundary& pdStandardFields::stdUBF()
{
    return stdU_.boundaryFieldRef();
}

//template <class pdParcel>
inline volScalarField::Boundary& pdStandardFields::transTBF()
{
    return transT_.boundaryFieldRef();
}*/


//template <class pdParcel>
// inline volScalarField::Boundary& pdStandardFields::boundaryRadiationTBF()
// {
//     return boundaryRadiationT_.boundaryFieldRef();
// }



//template <class pdParcel>
inline const volScalarField& pdStandardFields::q() const
{
    return q_;
}


//template <class pdParcel>
inline const volVectorField& pdStandardFields::fD() const
{
    return fD_;
}


//template <class pdParcel>
inline const volScalarField& pdStandardFields::rhoN() const
{
    return rhoN_;
}


//template <class pdParcel>
inline const volScalarField& pdStandardFields::rhoM() const
{
    return rhoM_;
}

//template <class pdParcel>
inline const volScalarField& pdStandardFields::pdRhoN() const
{
    return pdRhoN_;
}

//template <class pdParcel>
inline const volScalarField& pdStandardFields::linearKE() const
{
    return linearKE_;
}


//template <class pdParcel>
inline const volScalarField& pdStandardFields::rotationalE() const
{
    return rotationalE_;
}


//template <class pdParcel>
inline const volScalarField& pdStandardFields::rotationalDof() const
{
    return rotationalDof_;
}

inline const volScalarField& pdStandardFields::vibrationalE() const
{
    return vibrationalE_;
}


//template <class pdParcel>
inline const volScalarField& pdStandardFields::vibrationalDof() const
{
    return vibrationalDof_;
}

//template <class pdParcel>
inline const volVectorField& pdStandardFields::momentum() const
{
    return momentum_;
}
/*
//template <class pdParcel>
inline const volVectorField& pdStandardFields::stdU() const
{
    return stdU_;
}

//template <class pdParcel>
inline const volScalarField& pdStandardFields::transT() const
{
    return transT_;
}*/

// //template <class pdParcel>
// inline const volScalarField& pdStandardFields::boundaryRadiationT() const
// {
//     return boundaryRadiationT_;
// }



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
