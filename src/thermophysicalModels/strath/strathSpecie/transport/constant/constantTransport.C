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

\*---------------------------------------------------------------------------*/

#include "constantTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::constantTransport<Thermo>::constantTransport(Istream& is)
:
    Thermo(is),
    mu_(readScalar(is)),
    kappa_(readScalar(is)),
    kappave_(readScalar(is))
{
    is.check("constantTransport::constantTransport(Istream& is)");
}


template<class Thermo>
Foam::constantTransport<Thermo>::constantTransport(const dictionary& dict)
:
    Thermo(dict),
    mu_
    (
        dict.subDict("transport").subDict("constant")
            .lookupOrDefault<scalar>("mu", 0.0)
    ),
    kappa_
    (
        dict.subDict("transport").subDict("constant")
            .lookupOrDefault<scalar>("kappa", 0.0)
    ),
    kappave_
    (
        dict.subDict("transport").subDict("constant")
            .lookupOrDefault<scalar>("kappave", 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::constantTransport<Thermo>::constantTransport::write
(
    Ostream& os
) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dictTransport("transport");
    dictTransport.subDict("constant").add("mu", mu_);
    dictTransport.subDict("constant").add("kappa", kappa_);
    dictTransport.subDict("constant").add("kappave", kappave_);
    os  << indent << dictTransport.dictName() << dictTransport;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream&
Foam::operator<<(Ostream& os, const constantTransport<Thermo>& ct)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os << tab << ct.mu_ << tab << ct.kappa_ << tab << ct.kappave_;

    os.check("Ostream& operator<<(Ostream&, const constantTransport&)");

    return os;
}


// ************************************************************************* //
