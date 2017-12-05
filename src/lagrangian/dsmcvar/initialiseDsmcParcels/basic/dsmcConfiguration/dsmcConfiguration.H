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
    dsmcConfiguration

Description


SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef dsmcConfiguration_H
#define dsmcConfiguration_H

#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "Random.H"
#include "dsmcParcel.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcConfiguration Declaration
\*---------------------------------------------------------------------------*/

class dsmcConfiguration
{

protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        const dictionary& dsmcInitialiseDict_;

        Random& rndGen_;

        label nParcelsAdded_;

    // Private Member Functions

//         vector equipartitionLinearVelocity
//         (
//             scalar temperature,
//             scalar mass
//         );

public:

    //- Runtime type information
    TypeName("dsmcConfiguration");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcConfiguration,
            dictionary,
            (
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (cloud, dict)
        );
   
    // Constructors

        //- Construct from components
        dsmcConfiguration
        (
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcConfiguration> New
        (
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcConfiguration();


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        const label& nParcelsAdded() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
