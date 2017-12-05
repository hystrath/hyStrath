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
    polyField

Description

SourceFiles
    polyFieldI.H
    polyField.C
    polyFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyField_H
#define polyField_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMolecule.H"
#include "propertyField.H"
#include "writeTimeData.H"
#include "selectSiteIds.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyField Declaration
\*---------------------------------------------------------------------------*/

class polyField
{

protected:

    // Protected data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        Time& time_;

        fileName casePath_;

        fileName timePath_;

        bool measureInterForces_;
        bool measureInterForcesSites_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("polyField");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            polyField,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                polyMoleculeCloud& molCloud,
                const dictionary& dict
            ),
            (t, mesh, molCloud, dict)
        );
   
    // Constructors

        //- Construct from components
        polyField
        (
			Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<polyField> New
        (
			Time& t,
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const dictionary& dict
        );



    // Destructor

        virtual ~polyField();


    // Member Functions


        const fileName& casePath() const;
        fileName& casePath();

        const fileName& timePath() const;
        fileName& timePath();

        const bool& measureInterForces() const;
        bool& measureInterForces();

        const bool& measureInterForcesSites() const;
        bool& measureInterForcesSites();

        //- create initial field
        virtual void createField() = 0;

        //- calculate field
        virtual void calculateField() = 0;

        //- write field
        virtual void writeField() = 0;

        virtual void measureDuringForceComputation
        (
            polyMolecule* molI,
            polyMolecule* molJ
        ) = 0;

//         virtual void measureDuringForceComputation
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef
//         ) = 0;

        virtual void measureDuringForceComputationSite
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            label sI,
            label sJ
        ) = 0;

//         virtual void measureDuringForceComputationSite
//         (
//             polyMolecule* molReal,
//             polyReferredMolecule* molRef,
//             label sReal,
//             label sRef
//         ) = 0;


        virtual const propertyField& fields() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
