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
    polyDelFromZone

Description

SourceFiles
    polyDelFromZoneI.H
    polyDelFromZone.C
    polyDelFromZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyDelFromZone_H
#define polyDelFromZone_H

#include "polyMolsToDeleteModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyDelFromZone Declaration
\*---------------------------------------------------------------------------*/

class polyDelFromZone
:
	public polyMolsToDeleteModel
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;


        //- name of control zone
        word regionName_;
        label regionId_;

        bool oneSpecie_;

        label molId_;

    // Private Member Functions

        void findMolsToDel();

public:

    //- Runtime type information
    TypeName("polyDelFromZone");

    // Constructors

        //- Construct from components
        polyDelFromZone
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
		);

    // Destructor

 		~polyDelFromZone();


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
