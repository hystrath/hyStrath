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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const word& propertyField::fieldName() const
{
    return  fieldName_;
}

inline const volScalarField& propertyField::sField() const
{
    return  sField_;
}

inline  volScalarField& propertyField::sField() 
{
    return  sField_;
}

inline const volVectorField& propertyField::vField() const
{
    return  vField_;
}
inline  volVectorField& propertyField::vField() 
{
    return  vField_;
}
inline const volTensorField& propertyField::tField() const
{
    return  tField_;
}
inline  volTensorField& propertyField::tField() 
{
    return  tField_;
}

inline const scalar& propertyField::s() const
{
    return  s_;
}

inline  scalar& propertyField::s() 
{
    return  s_;
}

inline const vector& propertyField::v() const
{
    return  v_;
}
inline vector& propertyField::v() 
{
    return  v_;
}

inline const tensor& propertyField::t() const
{
    return  t_;
}

inline tensor& propertyField::t() 
{
    return  t_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
