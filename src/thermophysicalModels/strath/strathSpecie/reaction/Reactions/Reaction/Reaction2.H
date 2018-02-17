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

Class
    Foam::Reaction2

Description
    Simple extension of Reaction2Thermo to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    Reaction2I.H
    Reaction2.C

\*---------------------------------------------------------------------------*/

#ifndef Reaction2_H
#define Reaction2_H

#include "speciesTable.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

//#include "EnumString.h" // NEW VINCENT
//#include <cassert>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Reaction2Thermo>
class Reaction2;

template<class Reaction2Thermo>
inline Ostream& operator<<(Ostream&, const Reaction2<Reaction2Thermo>&);


// NEW VINCENT
enum ControllingTemperatureType {chargeExchange, dissociation, exchange, 
    impactDissociation, impactIonisation, associativeIonisation, 
    transrotational, vibrational};

/*---------------------------------------------------------------------------*\
                           Class Reaction2 Declaration
\*---------------------------------------------------------------------------*/

template<class Reaction2Thermo>
class Reaction2
:
    public Reaction2Thermo
{
protected:

    // Protected member functions

        //- Return string representation of the left of the reaction
        void reactionStrLeft(OStringStream& reaction) const;

        //- Return string representation of the right of the reaction
        void reactionStrRight(OStringStream& reaction) const;


public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedReactions;


    // Public data types

        //- Class to hold the specie index and its coefficients in the
        //  reaction rate expression
        struct specieCoeffs
        {
            label index;
            scalar stoichCoeff;
            scalar exponent;

            specieCoeffs()
            :
                index(-1),
                stoichCoeff(0),
                exponent(1)
            {}

            specieCoeffs(const speciesTable& species, Istream& is);

            bool operator==(const specieCoeffs& sc) const
            {
                return index == sc.index;
            }

            bool operator!=(const specieCoeffs& sc) const
            {
                return index != sc.index;
            }

            friend Ostream& operator<<(Ostream& os, const specieCoeffs& sc)
            {
                os  << sc.index << token::SPACE
                    << sc.stoichCoeff << token::SPACE
                    << sc.exponent;
                return os;
            }
        };


private:

    // Private data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;
        
        //- Controlling temperature
        ControllingTemperatureType controlT_; // NEW VINCENT


    // Private Member Functions

        //- Return string representation of reaction
        string reactionStr(OStringStream& reaction) const;

        //- Construct reaction thermo
        void setThermo(const HashPtrTable<Reaction2Thermo>& thermoDatabase);

        //- Disallow default bitwise assignment
        void operator=(const Reaction2<Reaction2Thermo>&);

        //- Return new reaction ID for un-named reactions
        label getNewReactionID();


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction2,
            Istream,
            (
                const speciesTable& species,
                const HashPtrTable<Reaction2Thermo>& thermoDatabase,
                Istream& is
            ),
            (species, thermoDatabase, is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction2,
            dictionary,
            (
                const speciesTable& species,
                const HashPtrTable<Reaction2Thermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );


    // Public classes

        //- Class used for the read-construction of PtrLists of reaction
        class iNew
        {
            const speciesTable& species_;
            const HashPtrTable<Reaction2Thermo>& thermoDatabase_;

        public:

            iNew
            (
                const speciesTable& species,
                const HashPtrTable<Reaction2Thermo>& thermoDatabase
            )
            :
                species_(species),
                thermoDatabase_(thermoDatabase)
            {}

            autoPtr<Reaction2> operator()(Istream& is) const
            {
                return autoPtr<Reaction2>
                (
                    Reaction2::New(species_, thermoDatabase_, is)
                );
            }
        };


    // Constructors

        //- Construct from components
        Reaction2
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs,
            const HashPtrTable<Reaction2Thermo>& thermoDatabase
        );

        //- Construct as copy given new speciesTable
        Reaction2(const Reaction2<Reaction2Thermo>&, const speciesTable& species);

        //- Construct from Istream
        Reaction2
        (
            const speciesTable& species,
            const HashPtrTable<Reaction2Thermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        Reaction2
        (
            const speciesTable& species,
            const HashPtrTable<Reaction2Thermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction2<Reaction2Thermo> > clone() const //
        {
            return autoPtr<Reaction2<Reaction2Thermo> >
            (
                new Reaction2<Reaction2Thermo>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction2<Reaction2Thermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction2<Reaction2Thermo> >
            (
                new Reaction2<Reaction2Thermo>(*this, species)
            );
        }


    // Selectors

        //- Return a pointer to new patchField created on freestore from input
        static autoPtr<Reaction2<Reaction2Thermo> > New
        (
            const speciesTable& species,
            const HashPtrTable<Reaction2Thermo>& thermoDatabase,
            Istream& is
        );

        //- Return a pointer to new patchField created on freestore from dict
        static autoPtr<Reaction2<Reaction2Thermo> > New
        (
            const speciesTable& species,
            const HashPtrTable<Reaction2Thermo>& thermoDatabase,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reaction2()
    {}


    // Member Functions

        // Access

            inline word& name();
            inline const word& name() const;

            // - Access to basic components of the reaction
            inline const List<specieCoeffs>& lhs() const;
            inline const List<specieCoeffs>& rhs() const;
            
            //- Controlling temperature of the reaction
            inline const ControllingTemperatureType& controlT() const;

            // - Access to gas components of the reaction
            virtual const List<specieCoeffs>& grhs() const;
            virtual const List<specieCoeffs>& glhs() const;

            // - Access to specie list
            const speciesTable& species() const;

            // - Access to gas specie list
            virtual const speciesTable& gasSpecies() const;

            //- Construct the left- and right-hand-side reaction coefficients
            void setLRhs
            (
                Istream&,
                const speciesTable&,
                List<specieCoeffs>& lhs,
                List<specieCoeffs>& rhs
            );


        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<< <Reaction2Thermo>
        (
            Ostream&,
            const Reaction2<Reaction2Thermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Reaction2I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Reaction2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
