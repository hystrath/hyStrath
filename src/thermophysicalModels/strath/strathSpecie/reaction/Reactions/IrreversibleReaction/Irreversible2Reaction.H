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
    Foam::Irreversible2Reaction

Description
    Simple extension of Reaction2 to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    Irreversible2Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Irreversible2Reaction_H
#define Irreversible2Reaction_H

#include "Reaction2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class Irreversible2Reaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class Irreversible2Reaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate k_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const Irreversible2Reaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&
        );


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        Irreversible2Reaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& reactionRate
        );

        //- Construct as copy given new speciesTable
        Irreversible2Reaction
        (
            const Irreversible2Reaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&,
            const speciesTable& species
        );

        //- Construct from Istream
        Irreversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        Irreversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction2<ReactionThermo> > clone() const
        {
            return autoPtr<Reaction2<ReactionThermo> >
            (
                new Irreversible2Reaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction2<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction2<ReactionThermo> >
            (
                new Irreversible2Reaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~Irreversible2Reaction()
    {}


    // Member Functions

        // Irreversible2Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Irreversible2Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
