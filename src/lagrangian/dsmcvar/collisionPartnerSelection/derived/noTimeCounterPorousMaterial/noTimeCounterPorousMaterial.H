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
 

Description

    no-time counter with reactions

SourceFiles
    noTimeCounterPorousMaterialI.H
    noTimeCounterPorousMaterial.C
    noTimeCounterPorousMaterialIO.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounterPorousMaterial_H
#define noTimeCounterPorousMaterial_H


#include "collisionPartnerSelection.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noTimeCounterPorousMaterial Declaration
\*---------------------------------------------------------------------------*/

class noTimeCounterPorousMaterial
:
    public collisionPartnerSelection
{

private:

    // Private data

//         dictionary propsDict_;




public:

    //- Runtime type information
    TypeName("noTimeCounterPorousMaterial");
    
    dictionary propsDict_;
    
    scalar porosity_;
    scalar temperature_;

    // Static data members

    // Constructors

        //- construct from components
        noTimeCounterPorousMaterial
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~noTimeCounterPorousMaterial();

    // Member Functions

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

        // Write





    // Member Operators

//         void operator=(const noTimeCounterPorousMaterial&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const noTimeCounterPorousMaterial&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
