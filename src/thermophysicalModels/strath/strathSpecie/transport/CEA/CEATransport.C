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
    temp_
    (
        dict.subDict("transport").isDict("CEA")
      ? dict.subDict("transport").subDict("CEA").lookup("temp") 
      : scalarList(4)
    ),
    mu_
    (
        dict.subDict("transport").isDict("CEA")
      ? dict.subDict("transport").subDict("CEA").lookup("visco")
      : List<CEATransportArray>(3)
    ),
    kappa_
    (
        dict.subDict("transport").isDict("CEA")
      ? dict.subDict("transport").subDict("CEA").lookup("kappa")
      : List<CEATransportArray>(3)
    ),
    eta_s_(dict.subDict("specie").lookupOrDefault<scalar>("eta_s", 1.2))
{
    if (not dict.subDict("transport").isDict("CEA"))
    {
        WarningInFunction
            << "Species: " << dict.dictName() << nl
            << "    transport/CEA subdictionary missing" << nl
            << "    CEA arrays temp, visco and kappa set to that of N2"
            << endl;
            
        temp_[0] = 200.0;
        temp_[1] = 1000.0;
        temp_[2] = 5000.0;
        temp_[3] = 15000.0;
        
        mu_[0][0] = 0.62526577;
        mu_[0][1] = -31.779652;
        mu_[0][2] = -1640.7983;
        mu_[0][3] = 1.7454992;
        
        mu_[1][0] = 0.87395209;
        mu_[1][1] = 561.52222;
        mu_[1][2] = -173948.09;
        mu_[1][3] = -0.39335958;
        
        mu_[2][0] = 0.88503551;
        mu_[2][1] = 909.02171;
        mu_[2][2] = -731290.61;
        mu_[2][3] = -0.53503838;
        
        kappa_[0][0] = 0.85439436;
        kappa_[0][1] = 105.73224;
        kappa_[0][2] = -12347.848;
        kappa_[0][3] = 0.47793128;
        
        kappa_[1][0] = 0.88407146;
        kappa_[1][1] = 133.57293;
        kappa_[1][2] = -11429.64;
        kappa_[1][3] = 0.24417019;
        
        kappa_[2][0] = 2.4176185;
        kappa_[2][1] = 8047.7749;
        kappa_[2][2] = 3105580.2;
        kappa_[2][3] = -14.517761;
    }

}


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
