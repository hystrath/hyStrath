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

#include "BlottnerEuckenTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport(Istream& is)
:
    Thermo(is),
    Ak_(readScalar(is)),
    Bk_(readScalar(is)),
    Ck_(readScalar(is))
{
    is.check("BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport(Istream&)");
}


template<class Thermo>
Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport(const dictionary& dict)
:
    Thermo(dict),
    Ak_(readScalar(dict.subDict("transport").subDict("BlottnerEucken").lookup("A"))),
    Bk_(readScalar(dict.subDict("transport").subDict("BlottnerEucken").lookup("B"))),
    Ck_(readScalar(dict.subDict("transport").subDict("BlottnerEucken").lookup("C")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::BlottnerEuckenTransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.subDict("BlottnerEucken").add("A", Ak_);
    dict.subDict("BlottnerEucken").add("B", Bk_);
    dict.subDict("BlottnerEucken").add("C", Ck_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    os << static_cast<const Thermo&>(bet) << tab << bet.Ak_ << tab << bet.Bk_ 
       << tab << bet.Ck_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const BlottnerEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
