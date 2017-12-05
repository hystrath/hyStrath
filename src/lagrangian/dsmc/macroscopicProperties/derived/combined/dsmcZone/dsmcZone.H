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
    dsmcZone

Description


SourceFiles
    dsmcZoneI.H
    dsmcZone.C
    dsmcZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcZone_H
#define dsmcZone_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
// #include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcZone Declaration
\*---------------------------------------------------------------------------*/

class dsmcZone
:
	public dsmcField
{

private:

    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

//         propertyField fields_;

//         autoPtr<binModel> binModel_;

        word regionName_;
        label regionId_;
        word fieldName_;
        scalar totalVolume_;
        
        labelList typeIds_;
        label timeIndex_;
        scalar averagingCounter_;        
        
        scalar mols_;
        scalar dsmcMols_;
        scalar molsInt_;
        scalar molsElec_;
        scalar mass_;
        scalar mcc_;
        vector mom_;
        vector UCollected_;
        scalar rotationalEMean_;
        scalar rotationalDofMean_;
        
        // pressure measurements
        
//         vectorField UMean_;
//         scalarField nD_;
//         vectorField Ucollected_;
//         scalarField massCollected_;
//         scalarField nParcels_;
        scalar muu_;
        scalar muv_;
        scalar muw_;
        scalar mvv_;
        scalar mvw_;
        scalar mww_;

        scalar mccu_;
        scalar mccv_;
        scalar mccw_;

        scalar eu_;
        scalar ev_;
        scalar ew_;
        scalar e_;
        
        label stepIndex_;
        scalarField speciesMols_;
        scalarField mccSpecies_;
        scalarField electronicETotal_;
        scalarField nParticlesGroundElectronicState_;
        scalarField nParticlesFirstElectronicState_;
        scalarField vDof_;
        scalarField mfp_;
        scalarField mcr_; // mean collision rate
        label nSteps_;
        
        scalarField N_;
        scalarField rhoN_;
        scalarField rhoM_;
        vectorField UMean_;
        vectorField UCAM_;        
        scalarField translationalTemperature_;
        scalarField rotationalTemperature_;
        scalarField vibrationalTemperature_;
        scalarField electronicTemperature_;
        scalarField overallTemperature_;
        scalarField scalarPressure_;
        tensorField pField_;
        tensorField tauField_;
        vectorField qField_;
        vectorField qInternalField_;
        vectorField qTranslationalField_;
        vectorField vibrationalModeTemperatures_;
        scalarField meanFreePath_;
        scalarField meanCollisionRate_;
        scalarField meanCollisionTime_;
        scalarField meanCollisionTimeTimeStepRatio_;
        scalarField Ma_;
        
        List<scalarField> vibrationalETotal_;
       
        List<bool> outputField_;

        bool instantaneous_;
        
        bool averagingAcrossManyRuns_;
        
        
    // Private Member Functions

        void readIn();

        void writeOut();

public:

    //- Runtime type information
    TypeName("dsmcZone");

    // Constructors

    //- Construct from components
    dsmcZone
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    // Destructor

    ~dsmcZone();


    // Member Functions

    //- create initial field
    void createField();

    //- calculate field
    void calculateField();

        //- return field
//         const volScalarField& densityField() const;
//         const propertyField& fields() const;
        //- write field
    void writeField();

    void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
