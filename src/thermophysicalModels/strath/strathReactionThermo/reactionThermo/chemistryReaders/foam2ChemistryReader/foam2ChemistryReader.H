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

Class
    Foam::foam2ChemistryReader

Description
    Chemistry reader for OpenFOAM format

SourceFiles
    foam2ChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef foam2ChemistryReader_H
#define foam2ChemistryReader_H

#include "chemistry2Reader.H"
#include "fileName.H"
#include "typeInfo.H"
#include "HashPtrTable.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class foam2ChemistryReader
:
    public chemistry2Reader<ThermoType>
{
        //- Chemistry/reactions dictionary
        dictionary chemDict_;

        //- Thermo properties dictionary
        dictionary thermoDict_;

        //- Table of species
        speciesTable& speciesTable_;
        
        //- Table of the thermodynamic data given in the foamChemistry file
        HashPtrTable<ThermoType> speciesThermo_;

        //- List of the reactions
        Reaction2List<ThermoType> reactions_;


    // Private Member Functions

        //- Set the species list
        speciesTable& setSpecies(const dictionary& dict, speciesTable& species);

        //- Disallow default bitwise copy construct
        foam2ChemistryReader(const foam2ChemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const foam2ChemistryReader&);
        
        //- Remove species present in the thermo dictionary and that are not
        //  listed in the table of species
        void removeAbsentSpecies();


public:

    //- Runtime type information
    TypeName("foam2ChemistryReader");


    // Constructors

        //- Construct from foamChemistry and thermodynamics file names
        foam2ChemistryReader
        (
            const fileName& reactionsFileName,
            speciesTable& species,
            const fileName& thermoFileName
        );

        //- Construct by getting the foamChemistry and thermodynamics file names
        //  from dictionary
        foam2ChemistryReader
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~foam2ChemistryReader()
    {}


    // Member functions

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of the thermodynamic data given in the foamChemistry file
        const HashPtrTable<ThermoType>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions
        const Reaction2List<ThermoType>& reactions() const
        {
            return reactions_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "foam2ChemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
