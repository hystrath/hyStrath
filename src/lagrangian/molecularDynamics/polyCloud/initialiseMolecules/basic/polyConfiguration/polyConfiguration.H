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
    polyConfiguration

Description


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef polyConfiguration_H
#define polyConfiguration_H

#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyConfiguration Declaration
\*---------------------------------------------------------------------------*/

class polyConfiguration
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        const dictionary& mdInitialiseDict_;

        label nMolsAdded_;

    // Private Member Functions
    
        // insert knowing from beforehand what cell info is
        void insertMolecule
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
        
        // insert and determine cell info in function
        void insertMolecule
        (
            const point& position,
            const label& id,         
            const bool& tethered,
            const bool& frozen,
            const scalar& temperature,
            const vector& bulkVelocity
        );
        

        void deleteMolecule
        (
            polyMolecule& mol
        );        
        
        vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

        vector equipartitionAngularMomentum
        (
            scalar temperature,
            label id
        );


public:

    //- Runtime type information
    TypeName("polyConfiguration");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyConfiguration,
            dictionary,
            (
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyConfiguration
        (
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyConfiguration> New
        (
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~polyConfiguration();


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        const label& mols() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
