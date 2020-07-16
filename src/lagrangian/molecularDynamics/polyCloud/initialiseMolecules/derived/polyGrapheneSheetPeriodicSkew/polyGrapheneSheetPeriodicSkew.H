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
    polyGrapheneSheetPeriodicSkew

Description


SourceFiles
    polyGrapheneSheetPeriodicSkewI.H
    polyGrapheneSheetPeriodicSkew.C
    polyGrapheneSheetPeriodicSkewIO.C


\*---------------------------------------------------------------------------*/

#ifndef polyGrapheneSheetPeriodicSkew_H
#define polyGrapheneSheetPeriodicSkew_H

#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyGrapheneSheetPeriodicSkew Declaration
\*---------------------------------------------------------------------------*/

class polyGrapheneSheetPeriodicSkew
:
	public polyConfiguration
{

private:

    // Private data

        vector breadthVector_;
        vector lengthVector_;

        bool initialiseVelocities_;

    // Private Member Functions



public:

    //- Runtime type information
    TypeName("polyGrapheneSheetPeriodicSkew");

    // Constructors

        //- Construct from components
        polyGrapheneSheetPeriodicSkew
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
//             const word& name
		);

    // Destructor

 		~polyGrapheneSheetPeriodicSkew();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
