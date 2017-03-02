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

#include "tabulated2KnTTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*template<class Thermo>
Foam::tabulated2KnTTransport<Thermo>::tabulated2KnTTransport(Istream& is)
:
    Thermo(is)
{    
    is.check("tabulated2KnTTransport<Thermo>::tabulated2KnTTransport(Istream&)");
    
    mu_ = interpolation2DTable<scalar>("constant/tabulatedKnTData/muTable");
    //mu_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
    kappa_ = interpolation2DTable<scalar>("constant/tabulatedKnTData/kappaTable");
    //kappa_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
}*/

template<class Thermo>
Foam::tabulated2KnTTransport<Thermo>::tabulated2KnTTransport(const dictionary& dict)
:
    Thermo(dict)
{
    word muFile_ = dict.subDict("transport").lookup("muFile");
    word kappaFile_ = dict.subDict("transport").lookup("kappaFile");
    mu_ = interpolation2DTable<scalar>("constant/tabulatedKnTData/" + muFile_);
    //mu_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
    kappa_ = interpolation2DTable<scalar>("constant/tabulatedKnTData/" + kappaFile_);
    //kappa_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::tabulated2KnTTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    /*dict.add("As", As_);
    dict.add("Ts", Ts_);*/
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const tabulated2KnTTransport<Thermo>& tt)
{
    operator<<(os, static_cast<const Thermo&>(tt));

    os.check("Ostream& operator<<(Ostream&, const tabulated2KnTTransport&)");

    return os;
}


// ************************************************************************* //
