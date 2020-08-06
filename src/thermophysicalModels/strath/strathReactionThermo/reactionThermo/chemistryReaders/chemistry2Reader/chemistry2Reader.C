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

#include "chemistry2Reader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistry2Reader<ThermoType> >
Foam::chemistry2Reader<ThermoType>::New
(
    const dictionary& thermoDict,
    speciesTable& species
)
{
    // Let the chemistry reader type default to CHEMKIN
    // for backward compatibility
    //word chemistry2ReaderTypeName("chemkinReader");
    word chemistry2ReaderTypeName("foam2ChemistryReader");

    // otherwise use the specified reader
    thermoDict.readIfPresent("chemistry2Reader", chemistry2ReaderTypeName);

    Info<< "Selecting chemistry2Reader " << chemistry2ReaderTypeName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistry2ReaderTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "chemistry2Reader::New(const dictionary&, speciesTable&)"
        )   << "Unknown chemistry2Reader type "
            << chemistry2ReaderTypeName << nl << nl
            << "Valid chemistry2Reader types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<chemistry2Reader<ThermoType> >
    (
        cstrIter()(thermoDict, species)
    );
}


// ************************************************************************* //
