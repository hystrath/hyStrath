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
    Foam::Reversible2Reaction

Description
    Simple extension of Reaction2 to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    Reversible2Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Reversible2Reaction_H
#define Reversible2Reaction_H

#include "Reaction2.H"

// NEW VINCENT 09/02/2016 *****************************************************
#include "DynamicList.H" 
#include "fundamentalConstants.H"
#include "physicoChemicalConstants.H"
// END NEW VINCENT 09/02/2016 *************************************************

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Reversible2Reaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class Reversible2Reaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate k_;
        
        // NEW VINCENT 09/02/2016 *********************************************
        //- Mixture number density table entries
        DynamicList<scalar> ni_; 
        
        //- Coefficients to determine Keq(nDMix): bounded interpolation
        DynamicList<scalar> A0_; 
        DynamicList<scalar> A1_; 
        DynamicList<scalar> A2_; 
        DynamicList<scalar> A3_; 
        DynamicList<scalar> A4_;
        // END NEW VINCENT 09/02/2016 *****************************************


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const Reversible2Reaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&
        );


public:

    //- Runtime type information
    TypeName("reversible");


    // Constructors

        //- Construct from components
        Reversible2Reaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& k,
            const DynamicList<scalar> ni, // NEW VINCENT 09/02/2016
            const DynamicList<scalar> A0,
            const DynamicList<scalar> A1,
            const DynamicList<scalar> A2,
            const DynamicList<scalar> A3,
            const DynamicList<scalar> A4
        );

        //- Construct as copy given new speciesTable
        Reversible2Reaction
        (
            const Reversible2Reaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&,
            const speciesTable& species
        );

        //- Construct from Istream
        Reversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        Reversible2Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<ReactionType<ReactionThermo> > clone() const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new Reversible2Reaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<ReactionType<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new Reversible2Reaction
                <ReactionType, ReactionThermo, ReactionRate>
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~Reversible2Reaction()
    {}


    // Member Functions

        // Reversible2Reaction rate coefficients

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
            //  Note this evaluates the forward rate constant and divides by
            //  the equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;
            
            // NEW VINCENT 09/02/2016 *****************************************
            //- Equilibrium constant calculation to evaluate the backward rate 
            //  constant
            virtual scalar Keq
            (
                const scalar p,
                const scalar T
            ) const;
            
            label findBounds
            (
                const scalar& xpt
            ) const
            {
                label bMax = ni_.size();
                label bCurrent = 0; 
                while (xpt > ni_[bCurrent+1] and bCurrent < bMax-2)
                {
                    bCurrent += 1;   
                }
                return bCurrent;
            }
            
            scalar boundedLinearInterpolation
            (
                const scalar& xpt,
                const scalar& x1, 
                const scalar& x2,
                const scalar& y1, 
                const scalar& y2
            ) const
            {
                if (xpt<=x1)
                {
                    return y1;
                }
                else if (xpt>=x2)
                {
                    return y2;
                }
                else
                {
                    return y1 + (xpt-x1)*(y2-y1)/(x2-x1);
                }
            }
            // END NEW VINCENT 09/02/2016 *************************************            


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Reversible2Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
