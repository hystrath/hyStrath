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
    Foam::chemistry2Reader

Description
    Abstract class for reading chemistry

SourceFiles
    chemistry2Reader.C

\*---------------------------------------------------------------------------*/

#ifndef chemistry2Reader_H
#define chemistry2Reader_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "Reaction2.H"
#include "Reaction2List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chemistry2Reader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistry2Reader
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistry2Reader(const chemistry2Reader&);

        //- Disallow default bitwise assignment
        void operator=(const chemistry2Reader&);


public:

    //- Runtime type information
    TypeName("chemistry2Reader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        chemistry2Reader()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistry2Reader,
            dictionary,
            (
                const dictionary& thermoDict,
                speciesTable& species
            ),
            (thermoDict, species)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<chemistry2Reader> New
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~chemistry2Reader()
    {}


    // Member Functions

        //- Return access to the list of species
        virtual const speciesTable& species() const = 0;

        //- Return access to the thermo packages
        virtual const HashPtrTable<ThermoType>& speciesThermo() const = 0;

        //- Return access to the list of reactions
        virtual const Reaction2List<ThermoType>& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistry2Reader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define make2ChemistryReader(MultiThermo)                                     \
    defineTemplateTypeNameAndDebug(chemistry2Reader<MultiThermo>, 0);         \
    defineTemplateRunTimeSelectionTable(chemistry2Reader<MultiThermo>, dictionary)


#define make2ChemistryReaderType(Reader, MultiThermo)                          \
    defineNamedTemplateTypeNameAndDebug(Reader<MultiThermo>, 0);               \
    chemistry2Reader<MultiThermo>::adddictionaryConstructorToTable<Reader<MultiThermo> > \
        add##Reader##MultiThermo##ConstructorToTable_

// for non-templated chemistry readers
#define add2ChemistryReaderType(Reader, MultiThermo)                          \
    defineTypeNameAndDebug(Reader, 0);                                        \
    chemistry2Reader<MultiThermo>::adddictionaryConstructorToTable<Reader>    \
        add##Reader##Thermo##ConstructorToTable_


// for templated chemistry readers
#define add2TemplateChemistryReaderType(Reader, MultiThermo)                  \
    defineNamedTemplateTypeNameAndDebug(Reader, 0);                           \
    chemistry2Reader<MultiThermo>::adddictionaryConstructorToTable<Reader>    \
        add##Reader##MultiThermo##ConstructorToTable_                         \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
