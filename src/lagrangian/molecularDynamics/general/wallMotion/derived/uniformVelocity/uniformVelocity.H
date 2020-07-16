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
    Wall is moved at constant speed.

Description

SourceFiles
    uniformVelocityI.H
    uniformVelocity.C
    uniformVelocityIO.C

\*---------------------------------------------------------------------------*/

#ifndef uniformVelocity_H
#define uniformVelocity_H


#include "wallMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniformVelocity Declaration
\*---------------------------------------------------------------------------*/

class uniformVelocity
:
    public wallMotion
{

private:

    // Private data

        dictionary propsDict_;

        vector velocity_;


public:

    //- Runtime type information
    TypeName("uniformVelocity");

    // Static data members

    // Constructors

        //- construct from components
        uniformVelocity
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~uniformVelocity();

    // Member Functions

//         void initialConfig
//         (
//             const scalar& density,
//             const radialDistribution& rdf,
//             const moleculeCloud& mD
//         );

        //- output data
//         void output(const Time& runTime, const polyMesh& mesh);


        // Access

        // Check

        // Edit

            void updateVelocity();

            const vector& velocity() const;

        // Write





    // Member Operators

//         void operator=(const uniformVelocity&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const uniformVelocity&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
