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
    Wall is started using a ramp increase in velocity.


Description

SourceFiles
    rampVelocityI.H
    rampVelocity.C
    rampVelocityIO.C

\*---------------------------------------------------------------------------*/

#ifndef rampVelocity_H
#define rampVelocity_H


#include "wallMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class rampVelocity Declaration
\*---------------------------------------------------------------------------*/

class rampVelocity
:
    public wallMotion
{

private:

    // Private data

        dictionary propsDict_;

        vector uInitial_;

        vector velocity_;

        vector uMax_;

        scalar tauT_;

        vector gradient_;

        scalar currentTimeElapsed_;

        scalar deltaTMD_;


public:

    //- Runtime type information
    TypeName("rampVelocity");

    // Static data members

    // Constructors

        //- construct from components
        rampVelocity
        (
            Time& time,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~rampVelocity();

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

//         void operator=(const rampVelocity&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const rampVelocity&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
