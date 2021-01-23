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

#include "SutherlandEuckenTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport(Istream& is)
:
    Thermo(is),
    As_(readScalar(is)),
    Ts_(readScalar(is)),
    eta_s_(readScalar(is))
{
    is.check
    (
        "SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport(Istream&)"
    );
}


template<class Thermo>
Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    As_(1.458e-6),
    Ts_(110.4),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{
    if (dict.subDict("transport").isDict("SutherlandEucken"))
    {
        As_ = dict.subDict("transport").subDict("SutherlandEucken")
            .lookupOrDefault<scalar>("As", 1.458e-6);
        Ts_ = dict.subDict("transport").subDict("SutherlandEucken")
            .lookupOrDefault<scalar>("Ts", 110.4);
    }
    else
    {
        WarningInFunction
            << "Species: " << dict.dictName() << nl
            << "    transport/SutherlandEucken subdictionary missing" << nl
            << "    Sutherland coefficients As and Ts set to that of Air" << nl
            << "    As = 1.458e-6, Ts = 110.4"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::SutherlandEuckenTransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dictTransport("transport");
    dictTransport.subDict("SutherlandEucken").add("As", As_);
    dictTransport.subDict("SutherlandEucken").add("Ts", Ts_);
    os  << indent << dictTransport.dictName() << dictTransport;

    dictionary dictSpecies("specie");
    dictSpecies.add("eta_s", eta_s_);
    os  << indent << dictSpecies.dictName() << dictSpecies;

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
    os << static_cast<const Thermo&>(st) << tab << st.As_ << tab << st.Ts_
       << tab << st.eta_s_;

    os.check
    (
        "Ostream& operator<<(Ostream&, "
            "const SutherlandEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
