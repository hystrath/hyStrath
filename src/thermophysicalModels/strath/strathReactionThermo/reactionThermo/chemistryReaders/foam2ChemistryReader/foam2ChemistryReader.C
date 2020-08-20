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

#include "foam2ChemistryReader.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::speciesTable& Foam::foam2ChemistryReader<ThermoType>::setSpecies
(
    const dictionary& dict,
    speciesTable& species
)
{
    wordList s(dict.lookup("species"));
    species.transfer(s);
    return species;
}


template<class ThermoType>
void Foam::foam2ChemistryReader<ThermoType>::removeAbsentSpecies()
{
    label j = 0;
    const label nSpecies = thermoDict_.keys().size();
    
    for(label i = 0; i<nSpecies; ++i)
    {
        bool speciesNotInMixture = true;
        
        forAll(speciesTable_, speciei)
        {
            if (thermoDict_.keys()[j] == speciesTable_[speciei])
            {
                speciesNotInMixture = false;
                break;
            }
        }
        
        if (speciesNotInMixture)
        {
            word speciesName = thermoDict_.keys()[j];
            thermoDict_.remove(speciesName);
            speciesThermo_.erase(speciesName);
        }
        else
        {
            j += 1;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foam2ChemistryReader<ThermoType>::foam2ChemistryReader
(
    const fileName& reactionsFileName,
    speciesTable& species,
    const fileName& thermoFileName
)
:
    chemistry2Reader<ThermoType>(),
    chemDict_
    (
        IFstream
        (
            fileName(reactionsFileName).expand()
        )()
    ),
    thermoDict_
    (
        IFstream
        (
            fileName(thermoFileName).expand()
        )()
    ),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    removeAbsentSpecies();
}


template<class ThermoType>
Foam::foam2ChemistryReader<ThermoType>::foam2ChemistryReader
(
    const dictionary& thermoDict,
    speciesTable& species
)
:
    chemistry2Reader<ThermoType>(),
    chemDict_
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryFile")).expand()
        )()
    ),
    thermoDict_
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
        )()
    ),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    removeAbsentSpecies();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
