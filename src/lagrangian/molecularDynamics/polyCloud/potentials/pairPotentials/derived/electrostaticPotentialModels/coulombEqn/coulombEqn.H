/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    coulombEqn


Description
    coulombEqn potential implemented in equation format

SourceFiles
    coulombEqnI.H
    coulombEqn.C
    coulombEqnIO.C

\*---------------------------------------------------------------------------*/

#ifndef coulombEqn_H
#define coulombEqn_H

#include "pairPotentialModel.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coulombEqn Declaration
\*---------------------------------------------------------------------------*/

class coulombEqn
:
    public pairPotentialModel
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;

        scalar constant_;

        // exclusion lists may go here

        
        
    // Private Member Functions



public:

    //- Runtime type information
    TypeName("coulombEqn");

    // Constructors

        //- Construct from components
        coulombEqn
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const reducedUnits& redUnits,
            const word& name,         
            const dictionary& dict
        );

    // Destructor

        ~coulombEqn();

    // Member Functions

        scalar unscaledEnergy(const scalar r) const;


        scalar force(const scalar r) const;
        
        scalar energy(const scalar r) const;        
        
//         void interaction (const scalar r, scalar& force, scalar& energy);
        //- Read dictionary
//         bool read
//         (
//             const dictionary& pairPotentialProperties,
//         );
        

        void write(const fileName& pathName);

        const dictionary& dict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
