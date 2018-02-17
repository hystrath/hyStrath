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

#include "multiThermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class MultiThermo, template<class> class Type>
const Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::tol_ = 1.0e-4; // Default OF: 1e-4; Gollan PhD09 suggests 1.0e-6

template<class MultiThermo, template<class> class Type>
const int Foam::species::multiThermo<MultiThermo, Type>::maxIter_ = 100;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
Foam::species::multiThermo<MultiThermo, Type>::multiThermo(Istream& is)
:
    MultiThermo(is)
{
    is.check("multiThermo<MultiThermo, Type>::multiThermo(Istream&)");
}


template<class MultiThermo, template<class> class Type>
Foam::species::multiThermo<MultiThermo, Type>::multiThermo(const dictionary& dict)
:
    MultiThermo(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
void Foam::species::multiThermo<MultiThermo, Type>::write(Ostream& os) const
{
    MultiThermo::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
Foam::Ostream& Foam::species::operator<<
(
    Ostream& os, const multiThermo<MultiThermo, Type>& st
)
{
    os  << static_cast<const MultiThermo&>(st);

    os.check("Ostream& operator<<(Ostream&, const multiThermo&)");
    return os;
}


// ************************************************************************* //
