/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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
        List<List<scalarField>> speciesRhoNIntBF_;
        List<List<scalarField>> speciesRhoNElecBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> speciesRhoMBF_;
        List<List<scalarField>> speciesLinearKEBF_;
        List<List<scalarField>> speciesMccBF_;
        List<List<vectorField>> speciesMomentumBF_;
        List<List<vectorField>> speciesUMeanBF_;
        List<List<scalarField>> speciesErotBF_;
        List<List<scalarField>> speciesZetaRotBF_;
        List<List<scalarField>> speciesEvibBF_;
        List<List<scalarField>> speciesEelecBF_;
        List<List<scalarField>> speciesqBF_;
        List<List<vectorField>> speciesfDBF_;

        //- species ID, vibrational energy mode, patch, face
        List<List<List<scalarField>>> speciesEvibModBF_;

        //- Number of stuck parcels on every boundary patch
        volScalarField::Boundary nParticlesOnStickingBoundaries_;

        //- Number of absorbed parcels on every boundary patch
        volScalarField::Boundary nAbsorbedParticles_;

        //- Boundary temperature field when dsmcFieldPatchBoundary is used
        volScalarField::Boundary boundaryT_;

        //- Boundary velocity field when dsmcFieldPatchBoundary is used
        volVectorField::Boundary boundaryU_;


    // Private Member Functions

        //- write the number of stuck parcels for sticking patches
        void writenStuckParticles();

        //- write the number of absorbed parcels for absorbing patches
        void writenAbsorbedParticles();

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
        void updatenStuckParticlesOnPatch
        (
            const label patchi,
            const scalarList& pnStuckParticles
        );

        //- update the number of absorbed parcels on patchi/facei
        void updatenAbsorbedParticlesOnPatch
        (
            const label patchi,
            const label facei,
            const scalar nAbsorbedParticles
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

        //- set the number of stuck parcels on patchi
        void setBoundarynStuckParticles
        (
            const label patchi,
            const scalarList& pnStuckParticles
        );

        //- set the number of absorbed parcels on patchi
        void setBoundarynAbsorbedParticles
        (
            const label patchi,
            const scalarList& pnAbsorbedParticles
        );

        void setInitialConfig();

        void clean();

        void reset();


    // Access

      inline scalar speciesRhoNIntBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesRhoNElecBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesRhoNBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesRhoMBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesLinearKEBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesMccBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline vector speciesMomentumBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline vector speciesUMeanBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesErotBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesZetaRotBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesEvibBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesEvibModBF
      (
          const label speciei,
          const label modi,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesEelecBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline scalar speciesqBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;
      
      inline vector speciesfDBF
      (
          const label speciei,
          const label patchi,
          const label facei
      ) const;

      //- Return const access to the number of stuck parcels on patchi
      const fvPatchScalarField& pnStuckParticles(const label patchi) const;

      //- Return const access to the number of absorbed parcels on patchi
      const fvPatchScalarField& pnAbsorbedParticles(const label patchi) const;

      //- Return access to the temperature boundary values on patchi
      const fvPatchScalarField& boundaryT(const label patchi) const;

      //- Return access to the velocity boundary values on patchi
      const fvPatchVectorField& boundaryU(const label patchi) const;


    // Edit

      inline scalar& speciesRhoNIntBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesRhoNElecBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesRhoNBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesRhoMBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesLinearKEBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesMccBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline vector& speciesMomentumBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline vector& speciesUMeanBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesErotBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesZetaRotBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesEvibBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesEvibModBF
      (
          const label speciei,
          const label modi,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesEelecBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline scalar& speciesqBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );
      
      inline vector& speciesfDBF
      (
          const label speciei,
          const label patchi,
          const label facei
      );

      //- Return access to the number of stuck parcels on patchi
      fvPatchScalarField& pnStuckParticles(const label patchi);

      //- Return access to the number of absorbed parcels on patchi
      fvPatchScalarField& pnAbsorbedParticles(const label patchi);

      //- Return access to the temperature boundary values on patchi
      fvPatchScalarField& boundaryT(const label patchi);

      //- Return access to the velocity boundary values on patchi
      fvPatchVectorField& boundaryU(const label patchi);


    // Write

      //- write boundary measurements fields at writeTime
      void outputResults();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
