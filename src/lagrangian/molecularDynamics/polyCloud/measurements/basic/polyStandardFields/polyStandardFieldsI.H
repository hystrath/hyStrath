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


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::qBF()
{
    return q_.boundaryField();
}


//template <class dsmcParcel>
inline volVectorField::GeometricBoundaryField& dsmcStandardFields::fDBF()
{
    return fD_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::rhoNBF()
{
    return rhoN_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::rhoMBF()
{
    return rhoM_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::linearKEBF()
{
    return linearKE_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::rotationalEBF()
{
    return rotationalE_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::rotationalDofBF()
{
    return rotationalDof_.boundaryField();
}

inline volScalarField::GeometricBoundaryField& dsmcStandardFields::vibrationalEBF()
{
    return vibrationalE_.boundaryField();
}


//template <class dsmcParcel>
inline volScalarField::GeometricBoundaryField& dsmcStandardFields::vibrationalDofBF()
{
    return vibrationalDof_.boundaryField();
}

//template <class dsmcParcel>
inline volVectorField::GeometricBoundaryField& dsmcStandardFields::momentumBF()
{
    return momentum_.boundaryField();
}


//template <class dsmcParcel>
// inline volScalarField::GeometricBoundaryField& dsmcStandardFields::boundaryRadiationTBF()
// {
//     return boundaryRadiationT_.boundaryField();
// }



//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::q() const
{
    return q_;
}


//template <class dsmcParcel>
inline const volVectorField& dsmcStandardFields::fD() const
{
    return fD_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::rhoN() const
{
    return rhoN_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::rhoM() const
{
    return rhoM_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::dsmcRhoN() const
{
    return dsmcRhoN_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::linearKE() const
{
    return linearKE_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::rotationalE() const
{
    return rotationalE_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::rotationalDof() const
{
    return rotationalDof_;
}

inline const volScalarField& dsmcStandardFields::vibrationalE() const
{
    return vibrationalE_;
}


//template <class dsmcParcel>
inline const volScalarField& dsmcStandardFields::vibrationalDof() const
{
    return vibrationalDof_;
}

//template <class dsmcParcel>
inline const volVectorField& dsmcStandardFields::momentum() const
{
    return momentum_;
}

// //template <class dsmcParcel>
// inline const volScalarField& dsmcStandardFields::boundaryRadiationT() const
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
