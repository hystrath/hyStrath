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
    boundaryMeasurements

Description

SourceFiles
    boundaryMeasurementsI.H
    boundaryMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryMeasurements_H
#define boundaryMeasurements_H

#include "vector.H"
#include "fvMesh.H"
#include "pdParcel.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boundaryMeasurements Declaration
\*---------------------------------------------------------------------------*/
// class pdCloud;


class boundaryMeasurements
{
    // Private data

        const fvMesh& mesh_;

        pdCloud& cloud_;

        labelList typeIds_;

        // Outer list is typeId, middle is the patch, inner is the local patch face
        List<scalarField> rhoNIntBF_;
        List< List < scalarField > > rhoNBF_;
        List< List < scalarField > > rhoMBF_;
        List< List < vectorField > > JpBF_;
        List< List < scalarField > > rhoQBF_;
        List< List < scalarField > > wallQBF_;
        List< List < scalarField > > linearKEBF_;
        List< List < vectorField > > momentumBF_;
        List< List < vectorField > > UMeanBF_;
        List< List < scalarField > > rotationalEBF_;
        List< List < scalarField > > rotationalDofBF_;
        List< List < scalarField > > vibrationalEBF_;
        List< List < scalarField > > qBF_;
        List< List < vectorField > > fDBF_;

    // Private Member Functions



public:

    friend class pdParcel;

    // Constructors

        //- Construct from mesh and cloud
        boundaryMeasurements
        (
            const polyMesh&,
            pdCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (pdFoam)
        boundaryMeasurements
        (
            const polyMesh&,
            pdCloud& cloud,
            const bool& dummy
        );

    // Destructor

        ~boundaryMeasurements();

    // Member Functions

        void updateFields
        (
            pdParcel& p
        );

        void clean();

        // Access
        inline const List<scalarField>& rhoNIntBF() const;
        inline const List< List < scalarField > >& rhoNBF() const;
        inline const List< List < scalarField > >& rhoMBF() const;
        inline const List< List < vectorField > >& JpBF() const;
        inline const List< List < scalarField > >& rhoQBF() const;
        inline const List< List < scalarField > >& wallQBF() const;
        inline const List< List < scalarField > >& linearKEBF() const;
        inline const List< List < vectorField > >& momentumBF() const;
        inline const List< List < vectorField > >& UMeanBF() const;
        inline const List< List < scalarField > >& rotationalEBF() const;
        inline const List< List < scalarField > >& rotationalDofBF() const;
        inline const List< List < scalarField > >& vibrationalEBF() const;
        inline const List< List < scalarField > >& qBF() const;
        inline const List< List < vectorField > >& fDBF() const;

        // Check

        // Edit

        inline List<scalarField>& rhoNIntBF();
        inline List< List < scalarField > >& rhoNBF();
        inline List< List < scalarField > >& rhoMBF();
        inline List< List < vectorField > >& JpBF();
        inline List< List < scalarField > >& rhoQBF();
        inline List< List < scalarField > >& wallQBF();
        inline List< List < scalarField > >& linearKEBF();
        inline List< List < vectorField > >& momentumBF();
        inline List< List < vectorField > >& UMeanBF();
        inline List< List < scalarField > >& rotationalEBF();
        inline List< List < scalarField > >& rotationalDofBF();
        inline List< List < scalarField > >& vibrationalEBF();
        inline List< List < scalarField > >& qBF();
        inline List< List < vectorField > >& fDBF();

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators


    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
