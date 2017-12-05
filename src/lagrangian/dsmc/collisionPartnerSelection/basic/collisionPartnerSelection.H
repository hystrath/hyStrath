/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    An abstract class for selecting collisions in a DSMC simulation

Description

SourceFiles
    collisionPartnerSelectionI.H
    collisionPartnerSelection.C
    collisionPartnerSelectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef collisionPartnerSelection_H
#define collisionPartnerSelection_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"
#include "fvMesh.H"
#include "dsmcParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class collisionPartnerSelection Declaration
\*---------------------------------------------------------------------------*/

class collisionPartnerSelection
{
    protected:

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        Random& rndGen_;


    private:


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("collisionPartnerSelection");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            collisionPartnerSelection,
            dictionary,
            (
                const polyMesh& mesh,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (mesh, cloud, dict)
        );

    // Static data members

    // Constructors

        collisionPartnerSelection
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Selectors

        static autoPtr<collisionPartnerSelection> New
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        virtual ~collisionPartnerSelection();

    // Member Functions

            virtual void initialConfiguration() = 0;

            virtual void collide() = 0;

        // Access

        // Check

        // Edit

        // Write

            //- output data

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "collisionPartnerSelectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
