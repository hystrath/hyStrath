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
#include "vectorList.H"
#include "fvMesh.H"
#include "dsmcParcel.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class boundaryMeasurements Declaration
\*---------------------------------------------------------------------------*/

//- Forward declaration of friend classes
class dsmcParcel;


class boundaryMeasurements
{

private:

    // Private data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        labelList typeIds_;
        
        //- Outer list is typeId, middle is the patch, inner is the local patch
        //  face
        List<List<scalarField>> rhoNIntBF_;
        List<List<scalarField>> rhoNElecBF_;
        List<List<scalarField>> rhoNBF_;
        List<List<scalarField>> rhoMBF_;
        List<List<scalarField>> linearKEBF_;
        List<List<scalarField>> mccSpeciesBF_;
        List<List<vectorField>> momentumBF_;
        List<List<vectorField>> UMeanBF_;
        List<List<scalarField>> rotationalEBF_;
        List<List<scalarField>> rotationalDofBF_;
        List<List<scalarField>> vibrationalEBF_;
        List<List<scalarField>> electronicEBF_;
        List<List<scalarField>> qBF_;
        List<List<vectorField>> fDBF_;
        
        //- Number of stuck parcels on every boundary patch
        volScalarField::Boundary nParcelsOnStickingBoundaries_;
        
        //- Number of absorbed parcels on every boundary patch
        volScalarField::Boundary nAbsorbedParcels_;
        
        //- Boundary temperature field when dsmcFieldPatchBoundary is used
        volScalarField::Boundary boundaryT_;
        
        //- Boundary velocity field when dsmcFieldPatchBoundary is used
        volVectorField::Boundary boundaryU_;


    // Private Member Functions
    
        //- write the number of stuck parcels for sticking patches
        void writenStuckParcels();
        
        //- write the number of absorbed parcels for absorbing patches
        void writenAbsorbedParcels();
        
        //- writePatchFields()
        void writePatchFields();
    


public:

    // Constructors

        //- Construct from mesh and cloud
        boundaryMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (dsmcFoam)
        boundaryMeasurements
        (
            const polyMesh&,
            dsmcCloud& cloud,
            const bool& dummy
        );


    // Destructor

        ~boundaryMeasurements();


    // Member Functions

        void updateFields(dsmcParcel& p);
        
        //- update the number of stuck parcels on patchi
        void updatenStuckParcelOnPatch
        (
            const label patchi,
            const scalarList& pnStuckParcels
        );
        
        //- update the number of absorbed parcels on patchi
        void updatenAbsorbedParcelOnPatch
        (
            const label patchi,
            const label facei
        );
        
        //- set the boundaryT field on patchi
        void setBoundaryT
        (
            const label patchi,
            const scalarList& pboundaryT
        );
        
        //- set the boundaryU field on patchi
        void setBoundaryU
        (
            const label patchi,
            const vectorList& pboundaryU
        );
        
        void clean();

        void reset();


    // Access
    
        const List<List<scalarField>>& rhoNIntBF() const;
        const List<List<scalarField>>& rhoNElecBF() const;
        const List<List<scalarField>>& rhoNBF() const;
        const List<List<scalarField>>& rhoMBF() const;
        const List<List<scalarField>>& linearKEBF() const;
        const List<List<scalarField>>& mccSpeciesBF() const;
        const List<List<vectorField>>& momentumBF() const;
        const List<List<vectorField>>& UMeanBF() const;
        const List<List<scalarField>>& rotationalEBF() const;
        const List<List<scalarField>>& rotationalDofBF() const;
        const List<List<scalarField>>& vibrationalEBF() const;
        const List<List<scalarField>>& electronicEBF() const;
        const List<List<scalarField>>& qBF() const;
        const List<List<vectorField>>& fDBF() const;
        
        //- Return const access to the number of stuck parcels on patchi
        const fvPatchScalarField& pnStuckParcels(const label& patchi) const;
        
        //- Return const access to the number of absorbed parcels on patchi
        const fvPatchScalarField& pnAbsorbedParcels(const label& patchi) const;
        
        //- Return access to the temperature boundary values on patchi
        const fvPatchScalarField& boundaryT(const label& patchi) const;
        
        //- Return access to the velocity boundary values on patchi
        const fvPatchVectorField& boundaryU(const label& patchi) const;


    // Check
    
    // Edit
    
        List<List<scalarField>>& rhoNIntBF();
        List<List<scalarField>>& rhoNElecBF();
        List<List<scalarField>>& rhoNBF();
        List<List<scalarField>>& rhoMBF();
        List<List<scalarField>>& linearKEBF();
        List<List<scalarField>>& mccSpeciesBF();
        List<List<vectorField>>& momentumBF();
        List<List<vectorField>>& UMeanBF();
        List<List<scalarField>>& rotationalEBF();
        List<List<scalarField>>& rotationalDofBF();
        List<List<scalarField>>& vibrationalEBF();
        List<List<scalarField>>& electronicEBF();
        List<List<scalarField>>& qBF();
        List<List<vectorField>>& fDBF();
        
        //- Return access to the number of stuck parcels on patchi
        fvPatchScalarField& pnStuckParcels(const label& patchi);
        
        //- Return access to the number of absorbed parcels on patchi
        fvPatchScalarField& pnAbsorbedParcels(const label& patchi);
        
        //- Return access to the temperature boundary values on patchi
        fvPatchScalarField& boundaryT(const label& patchi);
        
        //- Return access to the velocity boundary values on patchi
        fvPatchVectorField& boundaryU(const label& patchi);


    // Write
    
      //- write boundary measurements fields at writeTime
      void outputResults();


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
