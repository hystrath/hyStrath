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

#include "constantTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::constantTransport<Thermo>::constantTransport(Istream& is)
:
    Thermo(is),
    mu_(readScalar(is)),
    eta_s_(readScalar(is))
{
    is.check("constantTransport::constantTransport(Istream& is)");
}


template<class Thermo>
Foam::constantTransport<Thermo>::constantTransport(const dictionary& dict)
:
    Thermo(dict),
    mu_(readScalar(dict.subDict("transport").subDict("constant").lookup("mu"))),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::constantTransport<Thermo>::constantTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dictTransport("transport");
    dictTransport.subDict("constant").add("mu", mu_);
    os  << indent << dictTransport.dictName() << dictTransport;
    
    dictionary dictSpecies("specie");
    dictSpecies.add("eta_s", eta_s_);
    os  << indent << dictSpecies.dictName() << dictSpecies;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const constantTransport<Thermo>& ct)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os << tab << ct.mu_ << tab << ct.eta_s_;

    os.check("Ostream& operator<<(Ostream&, const constantTransport&)");

    return os;
}


// ************************************************************************* //
