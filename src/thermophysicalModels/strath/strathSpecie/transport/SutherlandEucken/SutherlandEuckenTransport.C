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

#include "SutherlandEuckenTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport(Istream& is)
:
    Thermo(is),
    As_(readScalar(is)),
    Ts_(readScalar(is))
{
    is.check("SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport(Istream&)");
}


template<class Thermo>
Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(readScalar(dict.subDict("transport").subDict("SutherlandEucken").lookup("As"))),
    Ts_(readScalar(dict.subDict("transport").subDict("SutherlandEucken").lookup("Ts")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::SutherlandEuckenTransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.subDict("SutherlandEucken").add("As", As_);
    dict.subDict("SutherlandEucken").add("Ts", Ts_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SutherlandEuckenTransport<Thermo>& st
)
{
    os << static_cast<const Thermo&>(st) << tab << st.As_ << tab << st.Ts_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const SutherlandEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
