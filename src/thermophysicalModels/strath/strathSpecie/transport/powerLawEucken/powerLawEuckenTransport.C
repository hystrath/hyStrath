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

#include "powerLawEuckenTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport(Istream& is)
:
    Thermo(is),
    dref_(readScalar(is)),
    omega_(readScalar(is)),
    eta_s_(readScalar(is))
{
    is.check("powerLawEuckenTransport<Thermo>::powerLawEuckenTransport(Istream&)");
}


template<class Thermo>
Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport(const dictionary& dict)
:
    Thermo(dict),
    dref_(readScalar(dict.subDict("specie").lookup("diameter"))),
    omega_(readScalar(dict.subDict("specie").lookup("omega"))),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::powerLawEuckenTransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("specie");
    dict.add("diameter", dref_);
    dict.add("omega", omega_);
    dict.add("eta_s", eta_s_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const powerLawEuckenTransport<Thermo>& pet
)
{
    os << static_cast<const Thermo&>(pet) << tab << pet.dref_ << tab 
       << pet.omega_ << tab << pet.eta_s_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const powerLawEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
