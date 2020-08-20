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

#include "BlottnerEuckenTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport(Istream& is)
:
    Thermo(is),
    Ak_(readScalar(is)),
    Bk_(readScalar(is)),
    Ck_(readScalar(is)),
    eta_s_(readScalar(is))
{
    is.check
    (
        "BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport(Istream&)"
    );
}


template<class Thermo>
Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    Ak_(2.68e-2),
    Bk_(3.18e-1),
    Ck_(-1.13e1),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{
    if (dict.subDict("transport").isDict("BlottnerEucken"))
    {
        Ak_ = dict.subDict("transport").subDict("BlottnerEucken")
            .lookupOrDefault<scalar>("A", 2.68e-2);
        Bk_ = dict.subDict("transport").subDict("BlottnerEucken")
            .lookupOrDefault<scalar>("B", 3.18e-1);
        Ck_ = dict.subDict("transport").subDict("BlottnerEucken")
            .lookupOrDefault<scalar>("C", -1.13e1);    
    }
    else
    {
        WarningInFunction
            << "Species: " << dict.dictName() << nl
            << "    transport/BlottnerEucken subdictionary missing" << nl
            << "    BlottnerEucken coefficients A, B and C set to that of N2"
            << nl << "    A = 2.68e-2, B = 3.18e-1, C = -1.13e1"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::BlottnerEuckenTransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dictTransport("transport");
    dictTransport.subDict("BlottnerEucken").add("A", Ak_);
    dictTransport.subDict("BlottnerEucken").add("B", Bk_);
    dictTransport.subDict("BlottnerEucken").add("C", Ck_);
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
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    os << static_cast<const Thermo&>(bet) << tab << bet.Ak_ << tab << bet.Bk_
       << tab << bet.Ck_ << tab << bet.eta_s_;

    os.check
    (
        "Ostream& operator<<(Ostream&, const BlottnerEuckenTransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
