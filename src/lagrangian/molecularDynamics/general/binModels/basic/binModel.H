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
    An abstract model for moving/controlling MD walls

Description

SourceFiles
    binModelI.H
    binModel.C
    binModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef binModel_H
#define binModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class binModel Declaration
\*---------------------------------------------------------------------------*/

class binModel
{

protected:


        const polyMesh& mesh_;


private:

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("binModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            binModel,
            dictionary,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Static data members

    // Constructors

        //- Construct from components
        binModel
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<binModel> New
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        virtual ~binModel();

    // Member Functions

            virtual label isPointWithinBin
            (
                const vector& rI,
                const label& cellI
            ) = 0;

            virtual scalarField binPositions() = 0;

            virtual vectorField bins() = 0;

            virtual const label& nBins() const = 0;

            virtual scalar binVolume(const label& n) = 0;

//             virtual scalar binWidth() = 0;

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

#include "binModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
