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

Class
    Foam::Reaction2List

Description
    List of templated reactions

SourceFiles
    Reaction2List.C

\*---------------------------------------------------------------------------*/

#ifndef Reaction2List_H
#define Reaction2List_H

#include "PtrList.H"
#include "SLPtrList.H"
#include "speciesTable.H"
#include "HashPtrTable.H"
#include "Reaction2.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class Reaction2List Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Reaction2List
:
    public SLPtrList<Reaction2<ThermoType> >
{
    // Private data

        //- Reference to the table of species
        const speciesTable& species_;

        //- Reference to the thermo database
        const HashPtrTable<ThermoType>& thermoDb_;

        //- The dictionary used for construction
        const dictionary dict_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Reaction2List&);


public:

    // Constructors

        //- Construct null
        Reaction2List
        (
            const speciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase
        );

        //- Construct from dictionary
        Reaction2List
        (
            const speciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct from file using (Istream)
        Reaction2List
        (
            const speciesTable& species,
            const HashPtrTable<ThermoType>& thermoDatabase,
            const fileName& fName
        );

        //- Construct copy
        Reaction2List(const Reaction2List<ThermoType>& reactions);


    //- Destructor
    ~Reaction2List();


    // Public Member Functions

        //- Read reactions from dictionary
        bool readReactionDict();

        //- Write
        void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Reaction2List.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
