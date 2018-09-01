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

#include "CEATransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::CEATransport<Thermo>::CEATransport(Istream& is)
:
    Thermo(is)
{
    forAll(temp_, i)
    {
        is >> temp_[i];
    }
    
    forAll(mu_, i)
    {
        forAll(mu_[i], j)
        {
            is >> mu_[i][j];
        }
    }
    
    forAll(kappa_, i)
    {
        forAll(kappa_[i], j)
        {
            is >> kappa_[i][j];
        }
    }
    
    is >> eta_s_;
    
    is.check("CEATransport<Thermo>::CEATransport(Istream&)");
}


template<class Thermo>
Foam::CEATransport<Thermo>::CEATransport(const dictionary& dict)
:
    Thermo(dict),
    temp_(dict.subDict("transport").subDict("CEA").lookup("temp")),
    mu_(dict.subDict("transport").subDict("CEA").lookup("visco")),
    kappa_(dict.subDict("transport").subDict("CEA").lookup("kappa")),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::CEATransport<Thermo>::write(Ostream& os) const
{
    os  << this->advancedSpecie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dictTransport("transport");
    dictTransport.subDict("CEA").add("temp", temp_);
    dictTransport.subDict("CEA").add("visco", mu_);
    dictTransport.subDict("CEA").add("kappa", kappa_);
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
    const CEATransport<Thermo>& ceat
)
{
    os << static_cast<const Thermo&>(ceat) << tab; 
    
    forAll(ceat.temp_, i)
    {
        os << ceat.temp_[i] << tab;
    }
    
    os << nl << "    ";
    
    forAll(ceat.mu_, i)
    {
        forAll(ceat.mu_[i], j)
        {
            os << ceat.mu_[i][j] << tab;
        }
    }
    
    os << nl << "    ";
    
    forAll(ceat.kappa_, i)
    {
        forAll(ceat.kappa_[i], j)
        {
            os << ceat.kappa_[i][j] << tab;
        }
    }
    
    os << ceat.eta_s_ << endl;
    
    os.check
    (
        "Ostream& operator<<(Ostream&, const CEATransport<Thermo>&)"
    );

    return os;
}


// ************************************************************************* //
