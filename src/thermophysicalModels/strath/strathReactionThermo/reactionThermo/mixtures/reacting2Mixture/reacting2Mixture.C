/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "reacting2Mixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reacting2Mixture<ThermoType>::reacting2Mixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    speciesTable(),
    autoPtr<chemistry2Reader<ThermoType> >
    (
        chemistry2Reader<ThermoType>::New(thermoDict, *this)
    ),
    multi2ComponentMixture<ThermoType>
    (
        thermoDict,
        *this,
        autoPtr<chemistry2Reader<ThermoType> >::operator()().speciesThermo(),
        mesh
    ),
    PtrList<Reaction2<ThermoType> >
    (
        autoPtr<chemistry2Reader<ThermoType> >::operator()().reactions()
    )
{
    autoPtr<chemistry2Reader<ThermoType> >::clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::reacting2Mixture<ThermoType>::read(const dictionary& thermoDict)
{}


// ************************************************************************* //
