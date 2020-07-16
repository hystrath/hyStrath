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
    polyGrapheneSheet

Description


SourceFiles
    polyGrapheneSheetI.H
    polyGrapheneSheet.C
    polyGrapheneSheetIO.C


\*---------------------------------------------------------------------------*/

#ifndef polyGrapheneSheet_H
#define polyGrapheneSheet_H

#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyGrapheneSheet Declaration
\*---------------------------------------------------------------------------*/

class polyGrapheneSheet
:
	public polyConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties

        vector breadthVector_;
        vector lengthVector_;

        bool initialiseVelocities_;

    // Private Member Functions






public:

    //- Runtime type information
    TypeName("polyGrapheneSheet");

    // Constructors

        //- Construct from components
        polyGrapheneSheet
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
//             const word& name
		);

    // Destructor

 		~polyGrapheneSheet();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
