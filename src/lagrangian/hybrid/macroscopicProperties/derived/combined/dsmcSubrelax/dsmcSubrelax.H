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
    dsmcSubrelax

Description
    Geometric field for temperature, density, velocity, pressure.



SourceFiles
    dsmcSubrelaxI.H
    dsmcSubrelax.C
    dsmcSubrelaxIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcSubrelax_H
#define dsmcSubrelax_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcSubrelax Declaration
\*---------------------------------------------------------------------------*/

class dsmcSubrelax
:
	public dsmcField
{

private:

    // Private data
        dictionary propsDict_;
        
        scalar mfpReferenceTemperature_;
        word fieldName_;
        label itBeforeRelaxing_;
        volScalarField dsmcRhoN_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField overallT_;
        volScalarField q_;
        volScalarField tau_;
        volScalarField meanFreePath_;
        volScalarField mfpCellRatio_;
        volScalarField meanCollisionRate_;
        volScalarField meanCollisionTime_;
        volScalarField meanCollisionTimeTimeStepRatio_;
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
        label nTimeSteps_;
        labelList typeIds_;
        scalarField rhoNMean_;
        scalarField rhoNMeanInt_;
        scalarField rhoMMean_;
        scalarField linearKEMean_;
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
        vectorField momentumMean_;
        List<labelList> boundaryCells_;
        List<scalarField> vibrationalETotal_;
        List<scalarField> nParcels_;
        List<scalarField> vibT_;

        // effective degrees of freedom 
        List<scalarField> vDof_;
        List<scalarField> mfp_; // mean free path
        List<scalarField> mcr_; // mean collision rate
        

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
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;
        
        List< List < scalarField > > vibrationalEBF_;
        List< List < scalarField > > speciesRhoNBF_;
        List< List < scalarField > > vibTBF_;
        List< List < scalarField > > vDofBF_;

        // starting fields
        scalarField rhoNMean0_;
        scalarField rhoNMeanInt0_;
        scalarField rhoMMean0_;
        scalarField linearKEMean0_;
        scalarField rotationalEMean0_;
        scalarField rotationalDofMean0_;
        scalarField muu0_;
        scalarField muv0_;
        scalarField muw0_;
        scalarField mvv0_;
        scalarField mvw0_;
        scalarField mww0_;
        scalarField mcc0_;
        scalarField mccu0_;
        scalarField mccv0_;
        scalarField mccw0_;
        scalarField eu0_;
        scalarField ev0_;
        scalarField ew0_;
        scalarField e0_;
        scalarField nClassI0_;
        scalarField nClassII0_;
        scalarField nClassIII0_;
        vectorField momentumMean0_;
        List<scalarField> vibrationalETotal0_;
        List<scalarField> nParcels0_;
        // starting boundary measures
        List<scalarField> rhoNBF0_;
        List<scalarField> rhoMBF0_;
        List<scalarField> linearKEBF0_;
        List<scalarField> rotationalEBF0_;
        List<scalarField> rotationalDofBF0_;
        List<scalarField> qBF0_;
        List<scalarField> speciesRhoNIntBF0_;
        List<vectorField> momentumBF0_;
        List<vectorField> fDBF0_;
        List< List < scalarField > > vibrationalEBF0_;
        List< List < scalarField > > speciesRhoNBF0_;

        // load initial variables
        const word pName_;
        const word UName_;
        const word TtName_;
        const word TrName_;
        const word mutrName_;
        const word kappatrName_;
        const List<word> rhoNNames_;
        const List<word> TvNames_;
        const List<word> kappavNames_;
        const scalar nEq_;
        const List<scalar> thetaV_;
        const List<scalar> mass_;
        const List<scalar> rotDOF_;

        bool averagingAcrossManyRuns_;
        bool measureClassifications_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool stopSamplingReset_;

    // Private Member Functions

        void readIn();

        void writeOut();

        void writeFiles();


public:
    //- Runtime type information
    TypeName("dsmcSubrelax");

    // Constructors

        //- Construct from components
    dsmcSubrelax
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~dsmcSubrelax();


    // Member Functions

    //- create initial field
    void createField();

    //- calculate field
    void calculateField();

    // For dsmcSubrelax
    void setReset
    (
        bool& stopReset
    );
    void setInitialFields
    (
        const List<word>& molecules
    );

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

