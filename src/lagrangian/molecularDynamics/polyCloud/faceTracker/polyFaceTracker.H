/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    polyFaceTracker

Description

    Class that handles tracking of molecules across faces and patches of the 
    domain 
    
    Refreshes data every time-step, from the cloud
    
SourceFiles
    polyFaceTrackerI.H
    polyFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef polyFaceTracker_H
#define polyFaceTracker_H

#include "vector.H"
#include "fvMesh.H"
#include "polyMolecule.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyFaceTracker Declaration
\*---------------------------------------------------------------------------*/
// class polyMoleculeCloud;


class polyFaceTracker
{
    // Private data

        const fvMesh& mesh_;

        polyMoleculeCloud& molCloud_;

        List<scalarField> molIdFlux_;
        List<scalarField> massIdFlux_;
        List<scalarField> absMomIdFlux_;
        List<vectorField> momIdFlux_;


    // Private Member Functions


public:

    friend class polyMolecule;

    // Constructors

        //- Construct from mesh and timeInfo
        polyFaceTracker
        (
            const polyMesh&,
            polyMoleculeCloud& molCloud
        );


    // Destructor

        ~polyFaceTracker();

    // Member Functions

//         void syncProcessorPatches();

//         void setIdFields
//         (
//             const label& nId
//         );

        void updateFields
        (
            polyMolecule& mol
//             const label& startingCell,
//             const label& crossedFace,
//             const scalar& mass,
//             const vector& U,
//             const vector& force,
//             const vector& mom,
//             const scalar& pE,
//             const scalar& kE,
//             const scalar& energy,
//             const scalar& volume
        );

        void clean();

        // Access


            inline const List<scalarField>& molIdFlux() const;

            inline const List<scalarField>& massIdFlux() const;

            inline const List<scalarField>& absMomIdFlux() const;

            inline const List<vectorField>& momIdFlux() const;

//             inline const surfaceScalarField& totMolFlux() const;
// 
//             inline const surfaceScalarField& molFlux() const;
// 
//             inline const surfaceScalarField& massFlux() const;
// 
//             inline const surfaceVectorField& momFlux() const;
// 
//             inline const surfaceScalarField& pEFlux() const;
// 
//             inline const surfaceScalarField& kEFlux() const;
// 
//             inline const surfaceScalarField& energyFlux() const;
// 
//             inline const surfaceScalarField& volumeFlux() const;
// 
//             inline const surfaceVectorField& forceFlux() const;

//             inline const List<scalarField>& idFlux() const;

        // Check

        // Edit

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators

//         inline friend bool operator==
//         (
//             const polyFaceTracker& a,
//             const polyFaceTracker& b
//         );
// 
//         inline friend bool operator!=
//         (
//             const polyFaceTracker& a,
//             const polyFaceTracker& b
//         );

    // IOstream Operators

//         friend Istream& operator>>(Istream&, polyFaceTracker&);
//         friend Ostream& operator<<(Ostream&, const polyFaceTracker&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyFaceTrackerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
