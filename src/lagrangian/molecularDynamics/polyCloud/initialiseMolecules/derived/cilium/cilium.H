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
    cilium

Description


SourceFiles
    ciliumI.H
    cilium.C
    ciliumIO.C


\*---------------------------------------------------------------------------*/

#ifndef cilium_H
#define cilium_H

#include "polyConfiguration.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "polyMoleculeCloud.H"
#include "distributePoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cilium Declaration
\*---------------------------------------------------------------------------*/

class cilium
:
	public polyConfiguration
{

private:

    // Private data
        
    
        DynamicList<label> trackingNumbers_;
//         DynamicList<vector> positions_;        

        DynamicList<label> trackingNumbersA_;
        DynamicList<label> trackingNumbersB_;

    // Private Member Functions

        
        void insertMoleculeLocal
        (
            const point& position,
            const label cell,
            const label tetFace,
            const label tetPt, 
            const label& id, 
            const bool& tethered,
            const bool& frozen,
            const scalar& temperature,
            const vector& bulkVelocity
        );

        
public:

    //- Runtime type information
    TypeName("cilium");

    // Constructors

        //- Construct from components
        cilium
        (
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );

    // Destructor

        ~cilium();


    // Member Functions

        void setInitialConfiguration();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
