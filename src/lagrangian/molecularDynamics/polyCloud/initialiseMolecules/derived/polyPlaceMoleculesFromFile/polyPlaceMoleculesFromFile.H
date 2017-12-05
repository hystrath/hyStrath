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
    polyPlaceMoleculesFromFile

Description


SourceFiles
    polyPlaceMoleculesFromFileI.H
    polyPlaceMoleculesFromFile.C
    polyPlaceMoleculesFromFileIO.C


\*---------------------------------------------------------------------------*/

#ifndef polyPlaceMoleculesFromFile_H
#define polyPlaceMoleculesFromFile_H

#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyPlaceMoleculesFromFile Declaration
\*---------------------------------------------------------------------------*/

class polyPlaceMoleculesFromFile
:
	public polyConfiguration
{

private:

    // Private data

        //- subDictionary containing the properties
//         dictionary propsDict_;
//         fileName timePath_;
        
        word option_;

    // Private Member Functions


        void insertMolecule
        (
            const point& position,
            const label& cell,
            const label& id,
            const label& tetFace,
            const label& tetPt,         
            const bool& tethered,
            const bool& frozen,
            const scalar& phi,
            const scalar& theta,
            const scalar& psi,
//             const scalar& temperature,
            const vector& velocity
        );

        void fixedPropertiesFromDict();
        
        void fixedPropertiesFromFile();
        
        void multiIdFixedProperties();
        
        
public:

    //- Runtime type information
    TypeName("polyPlaceMoleculesFromFile");

    // Constructors

        //- Construct from components
        polyPlaceMoleculesFromFile
		(
            polyMoleculeCloud& molCloud,
            const dictionary& dict
//             const word& name
		);

    // Destructor

 		~polyPlaceMoleculesFromFile();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
