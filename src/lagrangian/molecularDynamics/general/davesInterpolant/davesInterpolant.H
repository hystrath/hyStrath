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


Description


SourceFiles
    davesInterpolantI.H
    davesInterpolant.C
    davesInterpolantIO.C

\*---------------------------------------------------------------------------*/

#ifndef davesInterpolant_H
#define davesInterpolant_H

#include "scalarField.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class davesInterpolant Declaration
\*---------------------------------------------------------------------------*/

class davesInterpolant
{
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        label nDataPoints_;


        scalarField density_;

        List< List<scalar> > strain_;
        List< List<scalar> > viscosity_;
        List< List<scalar> > slipLength_;




    // Private Member Functions

//         void setRadius();

public:

    // Static data members

    // Constructors

        //- Construct from dictionary
        davesInterpolant(const IOdictionary& dict);


    // Selectors

    // Destructor

        ~davesInterpolant();

    // Member Functions



    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "davesInterpolantI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
