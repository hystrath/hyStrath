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
    dsmcVolFields

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    dsmcVolFieldsI.H
    dsmcVolFields.C
    dsmcVolFieldsIO.C

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

private:

    // Private data
    
        dictionary propsDict_;
        
        label sampleInterval_;
        label sampleCounter_;
        scalar mfpReferenceTemperature_;
        word fieldName_;
        volScalarField dsmcRhoN_;
        volScalarField dsmcRhoNMean_;
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
        scalarField rhoNMean_;
        scalarField rhoNInstantaneous_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoNMeanInt_;
        scalarField molsElec_;
        scalarField rhoMMean_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMean_;
        scalarField linearKEMeanXnParticle_;
        scalarField rotationalEMean_;
        scalarField rotationalDofMean_;
        scalarField muu_;
        scalarField muv_;
        scalarField muw_;
        scalarField mvv_;
        scalarField mvw_;
        scalarField mww_;
        scalarField mcc_;
        scalarField mccu_;
        scalarField mccv_;
        scalarField mccw_;
        scalarField eu_;
        scalarField ev_;
        scalarField ew_;
        scalarField e_;
        scalarField totalvDof_;
        scalarField nClassI_;
        scalarField nClassII_;
        scalarField nClassIII_;
        scalarField collisionSeparation_;
        scalarField nColls_;
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;
        List<labelList> boundaryCells_;
        List< List < scalarField > > vibrationalETotal_;
        List<scalarField> electronicETotal_;
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;
        List<scalarField> vibT_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;

        //- Effective degrees of freedom 
        List<scalarField> vDof_;
        
        //- Mean free path
        List<scalarField> mfp_; 
        
        //- Mean collision rate
        List<scalarField> mcr_;
        
        //- Measured local collision rate 
        volScalarField cr_; 
        

        // boundary measurements
        
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
        
    // Private Member Functions

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
      void createField();

      //- calculate field
      void calculateField();
      
      //- reset field if mesh has changed
      void resetField();

      //- write field
      void writeField();

  //         const propertyField& fields() const;
      void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

