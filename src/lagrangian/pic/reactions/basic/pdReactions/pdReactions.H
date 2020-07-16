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
    pdReactions

Description

	Stores all the information for the pdReactions models defined within
    the pdReactionsDict, and selects & builds the models automatically.



\*---------------------------------------------------------------------------*/

#ifndef pdReactions_H
#define pdReactions_H


#include "List.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "polyMesh.H"
// #include "timeData.H"

#include "pdReaction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdReactions Declaration
\*---------------------------------------------------------------------------*/

class pdReactions
{

    // Private data

        //- The entire dictionary (containing multiple subDictionaries)

        Time& time_;

        IOdictionary chemReactDict_;

        // number of pdReactions

        label nReactions_;

        PtrList<entry> reactionsList_;
        List<word> reactionNames_;
        List<label> reactionIds_;
        List< autoPtr<pdReaction> > reactions_;

        // given two pair typeIds, one is returned the modelId number in the reactions_ data member
        List< List<label> > pairAddressing_;

        label counter_;

public:

    // Constructors

        //- Null Constructor
        pdReactions
        (
            Time& t,
            const polyMesh& mesh
        );

        //- Constructor for mdFOAM
        pdReactions
        (
            Time& t,
            const polyMesh& mesh,
            pdCloud& cloud
        );


        ~pdReactions();

    // Member Functions

        //- initial configuration
        void initialConfiguration();

        void outputData();

        // Access

            inline label returnModelId(pdParcel& p, pdParcel& q);

            inline List< autoPtr<pdReaction> >& reactions();
            inline const List< autoPtr<pdReaction> >& reactions() const;

            inline const List< List<label> >& pairModelAddressing() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdReactionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
