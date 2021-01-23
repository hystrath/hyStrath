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

#include "powerLawEuckenTransport.H"
#include "IOstreams.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport(Istream& is)
:
    Thermo(is),
    Tref_(273.0),
    dref_(readScalar(is)),
    omega_(readScalar(is)),
    muref_
    (
        15.0
      * sqrt
        (
            constant::mathematical::pi*constant::physicoChemical::k.value()
          * (this->W()*1.0e-3/constant::physicoChemical::NA.value())*Tref_
        )
       /
        (
            constant::mathematical::twoPi*sqr(dref_)
          * (5.0-2.0*omega_)*(7.0-2.0*omega_)
        )
    ),
    eta_s_(readScalar(is))
{
    is.check
    (
        "powerLawEuckenTransport<Thermo>::powerLawEuckenTransport(Istream&)"
    );
}


template<class Thermo>
Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    Tref_(273.0),
    dref_(readScalar(dict.subDict("specie").lookup("diameter"))),
    omega_(readScalar(dict.subDict("specie").lookup("omega"))),
    muref_
    (
        15.0
      * sqrt
        (
            constant::mathematical::pi*constant::physicoChemical::k.value()
          * (this->W()*1.0e-3/constant::physicoChemical::NA.value())*Tref_
        )
       /
        (
            constant::mathematical::twoPi*sqr(dref_)
          * (5.0-2.0*omega_)*(7.0-2.0*omega_)
        )
    ),
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
    os << static_cast<const Thermo&>(pet) << tab << pet.Tref_ << tab
       << pet.dref_ << tab << pet.omega_ << tab << pet.muref_ << tab
       << pet.eta_s_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const powerLawEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
