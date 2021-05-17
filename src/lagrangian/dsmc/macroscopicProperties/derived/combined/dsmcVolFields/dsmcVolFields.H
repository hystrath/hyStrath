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
    dsmcVolFields

Description
    Geometric field for temperature, density, velocity, pressure.

SourceFiles
    dsmcVolFields.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVolFields_H
#define dsmcVolFields_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dsmcVolFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcVolFields
:
	public dsmcField
{

protected:

    // Protected data

        dictionary propsDict_;

        label sampleInterval_;
        label sampleCounter_;
        scalar nTimeSteps_;
        
        scalar mfpTref_;
        word fieldName_;
        
        labelList speciesIds_;
        
        volScalarField dsmcN_;
        volScalarField dsmcNMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField Ttra_;
        volScalarField Trot_;
        volScalarField Tvib_;
        volScalarField Telec_;
        volScalarField Tov_;
        volScalarField q_;
        volScalarField tau_;
        volScalarField mfp_;
        volScalarField mfpToDx_;
        volScalarField DxToMfp_;
        volScalarField meanCollisionRate_;
        volScalarField measuredCollisionRate_;
        volScalarField meanCollisionTime_;
        volScalarField mctToDt_;
        volScalarField meanCollisionSeparation_;
        volScalarField SOF_;
        volScalarField Ma_;
        volScalarField classIDistribution_;
        volScalarField classIIDistribution_;
        volScalarField classIIIDistribution_;
        volScalarField densityError_;
        volScalarField velocityError_;
        volScalarField temperatureError_;
        volScalarField pressureError_;
        volVectorField UMean_;
        volVectorField fD_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;

        scalarField dsmcNCum_;
        scalarField nCum_;
        scalarField dsmcNElecLvlCum_;
        scalarField dsmcMCum_;
        scalarField mCum_;
        scalarField dsmcLinearKECum_;
        scalarField linearKECum_;
        scalarField dsmcErotCum_;
        scalarField dsmcZetaRotCum_;
        scalarField dsmcMuuCum_;
        scalarField dsmcMuvCum_;
        scalarField dsmcMuwCum_;
        scalarField dsmcMvvCum_;
        scalarField dsmcMvwCum_;
        scalarField dsmcMwwCum_;
        scalarField dsmcMccCum_;
        scalarField dsmcMccuCum_;
        scalarField dsmcMccvCum_;
        scalarField dsmcMccwCum_;
        scalarField dsmcEuCum_;
        scalarField dsmcEvCum_;
        scalarField dsmcEwCum_;
        scalarField dsmcECum_;
        scalarField zetaVib_;
        scalarField dsmcNClassICum_;
        scalarField dsmcNClassIICum_;
        scalarField dsmcNClassIIICum_;
        scalarField collisionSeparation_;
        scalarField dsmcNCollsCum_;
        vectorField dsmcMomentumCum_;
        vectorField momentumCum_;
        
        List<labelList> boundaryCells_;
        List<List<scalarField>> dsmcSpeciesEvibModCum_;
        List<scalarField> dsmcSpeciesEelecCum_;
        List<scalarField> dsmcNSpeciesCum_;
        List<scalarField> nSpeciesCum_;
        List<scalarField> dsmcMccSpeciesCum_;
        List<scalarField> speciesTvib_;
        List<scalarField> dsmcNGrndElecLvlSpeciesCum_;
        List<scalarField> dsmcN1stElecLvlSpeciesCum_;
        List<scalarField> speciesMfp_;
        List<scalarField> speciesMcr_;


      //- Boundary measurements
      
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> ErotBF_;
        List<scalarField> zetaRotBF_;
        List<scalarField> qBF_;
        List<scalarField> zetaVibBF_;
        List<scalarField> rhoNIntBF_;
        List<scalarField> rhoNElecBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;

        List<List<scalarField>> speciesEvibBF_;
        List<List<scalarField>> speciesEelecBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> speciesMccBF_;
        List<List<scalarField>> speciesTvibBF_;
        List<List<scalarField>> speciesZetaVibBF_;
        List<List<List<scalarField>>> speciesEvibModBF_;

        //- Wall unit vectors
        //  temporary storage, should be done in boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        bool averagingAcrossManyRuns_;
        bool measureClassifications_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool densityOnly_;
        bool measureHeatFluxShearStress_;

    // Protected Member Functions

        void readIn();

        void writeOut();

        void calculateWallUnitVectors();


public:

    //- Runtime type information
    TypeName("dsmcVolFields");

    // Constructors

        //- Construct from components
        dsmcVolFields
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

    ~dsmcVolFields();


    // Member Functions

      //- create initial field
      virtual void createField();

      //- calculate field
      virtual void calculateField();

      //- reset field if mesh has changed
      virtual void resetField();

      //- write field
      virtual void writeField();

      virtual void updateProperties(const dictionary& dict);

    // Access

      virtual scalar translationalT(const label cellI)
      {
          return Ttra_[cellI];
      }

      virtual scalar overallT(const label cellI)
      {
          return Tov_[cellI];
      }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

