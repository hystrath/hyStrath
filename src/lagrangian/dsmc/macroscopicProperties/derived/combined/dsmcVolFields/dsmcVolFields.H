/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
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
        scalar mfpReferenceTemperature_;
        word fieldName_;
        volScalarField dsmcN_;
        volScalarField dsmcNMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField electronicT_;
        volScalarField overallT_;
        volScalarField q_;
        volScalarField tau_;
        volScalarField meanFreePath_;
        volScalarField mfpCellRatio_;
        volScalarField cellMfpRatio_;
        volScalarField meanCollisionRate_;
        volScalarField meanCollisionTime_;
        volScalarField meanCollisionTimeTimeStepRatio_;
        //- Measured local collision rate
        volScalarField cr_;
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
        scalar nTimeSteps_;
        labelList typeIds_;
        scalarField dsmcNCum_;
        scalarField dsmcNInstantaneous_;
        scalarField nCum_;
        scalarField dsmcNWithRotDofCum_;
        scalarField dsmcNEleCum_;
        scalarField dsmcMCum_;
        scalarField mCum_;
        scalarField dsmcLinKinEnCum_;
        scalarField linKinEnCum_;
        scalarField dsmcRotEnCum_;
        scalarField dsmcRotDofCum_;
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
        scalarField totalvDof_;
        scalarField dsmcNClassICum_;
        scalarField dsmcNClassIICum_;
        scalarField dsmcNClassIIICum_;
        scalarField collisionSeparation_;
        scalarField dsmcNCollsCum_;
        vectorField dsmcMomentumCum_;
        vectorField momentumCum_;
        List<labelList> boundaryCells_;
        List< List < scalarField > > dsmcVibEnSpeciesModeCum_;
        List<scalarField> dsmcEleEnSpeciesCum_;
        List<scalarField> dsmcNSpeciesCum_;
        List<scalarField> nSpeciesCum_;
        List<scalarField> dsmcMccSpeciesCum_;
        List<scalarField> vibT_;
        List<scalarField> dsmcNGroundEleLvlSpeciesCum_;
        List<scalarField> dsmcNFirstEleLvlSpeciesCum_;

        //- Effective degrees of freedom
        List<scalarField> vDof_;

        //- Mean free path
        List<scalarField> mfp_;

        //- Mean collision rate
        List<scalarField> mcr_;




        //- Boundary measurements
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> rotationalEBF_;
        List<scalarField> rotationalDofBF_;
        List<scalarField> qBF_;
        List<scalarField> vibTxvDofBF_;
        List<scalarField> totalvDofBF_;
        List<scalarField> speciesRhoNIntBF_;
        List<scalarField> speciesRhoNElecBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;

        List<List<scalarField>> vibrationalEBF_;
        List<List<scalarField>> electronicEBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> mccSpeciesBF_;
        List<List<scalarField>> vibTBF_;
        List<List<scalarField>> vDofBF_;
        List<List<List<scalarField>>> evmsBF_;

        //- Wall unit vectors - temporary storage, should be done in boundaries/patch
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
          return translationalT_[cellI];
      }

      virtual scalar overallT(const label cellI)
      {
          return overallT_[cellI];
      }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

