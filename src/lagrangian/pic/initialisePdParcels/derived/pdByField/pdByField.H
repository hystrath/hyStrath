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
    pdByField

Description


SourceFiles
    pdByFieldI.H
    pdByField.C
    pdByFieldIO.C


\*---------------------------------------------------------------------------*/

#ifndef pdByField_H
#define pdByField_H

#include "pdAllConfigurations.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "pdCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pdByField Declaration
\*---------------------------------------------------------------------------*/

class pdByField
:
	public pdConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;
    List<scalarField>   numberDensities_;
    List<vectorField>   particleVelocity_;

    // Private Member Functions




public:

    //- Runtime type information
    TypeName("pdByField");

    // Constructors

        //- Construct from components
        pdByField
		(
            pdCloud& cloud,
            const dictionary& dict
		);

    // Destructor

 		~pdByField();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
