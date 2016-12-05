/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "perfect2Gas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2Gas<Specie>::perfect2Gas(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2Gas<Specie>::perfect2Gas
(
    const word& name,
    const perfect2Gas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::perfect2Gas<Specie> >
Foam::perfect2Gas<Specie>::clone() const
{
    return autoPtr<perfect2Gas<Specie> >(new perfect2Gas<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::perfect2Gas<Specie> >
Foam::perfect2Gas<Specie>::New(Istream& is)
{
    return autoPtr<perfect2Gas<Specie> >(new perfect2Gas<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::perfect2Gas<Specie> > Foam::perfect2Gas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<perfect2Gas<Specie> >(new perfect2Gas<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::perfect2Gas<Specie>::rho(scalar p, scalar T) const
{
    //Info << "rho: "<< p/(this->R()*T) << endl;
    return p/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfect2Gas<Specie>::psi(scalar, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfect2Gas<Specie>::Z(scalar, scalar) const
{
    return 1.0;
}


template<class Specie>
inline Foam::scalar Foam::perfect2Gas<Specie>::cpMcv(scalar, scalar) const
{
    return this->RR;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::perfect2Gas<Specie>::operator+=(const perfect2Gas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::perfect2Gas<Specie>::operator-=(const perfect2Gas<Specie>& pg)
{
    Specie::operator-=(pg);
}


template<class Specie>
inline void Foam::perfect2Gas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2Gas<Specie> Foam::operator+
(
    const perfect2Gas<Specie>& pg1,
    const perfect2Gas<Specie>& pg2
)
{
    return perfect2Gas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::perfect2Gas<Specie> Foam::operator-
(
    const perfect2Gas<Specie>& pg1,
    const perfect2Gas<Specie>& pg2
)
{
    return perfect2Gas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::perfect2Gas<Specie> Foam::operator*
(
    const scalar s,
    const perfect2Gas<Specie>& pg
)
{
    return perfect2Gas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::perfect2Gas<Specie> Foam::operator==
(
    const perfect2Gas<Specie>& pg1,
    const perfect2Gas<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
