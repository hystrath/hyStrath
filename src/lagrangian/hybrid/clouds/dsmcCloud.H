/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::dsmcCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    dsmcCloudI.H
    dsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCloud_H
#define dsmcCloud_H

#include "Cloud.H"

#include "dsmcParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "BinaryCollisionModel.H"

#include "dsmcBoundaries.H"
#include "dsmcAllConfigurations.H"
#include "dsmcFieldProperties.H"
#include "dsmcFaceTracker.H"
#include "dsmcControllers.H"
#include "dsmcStandardFields.H"

#include "dsmcReactions.H"
#include "collisionPartnerSelection.H"
// #include "reactInfoList.H"
#include "boundaryMeasurements.H"

#include "molsToDeleteHybrid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dsmcCloud Declaration
\*---------------------------------------------------------------------------*/


class dsmcCloud
:
    public Cloud<dsmcParcel>
{

private:

    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudName_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of unique instances of molecule types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<dsmcParcel*> > cellOccupancy_;

        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<dsmcParcel::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- controllers
        dsmcControllers controllers_;

        //- measurement fields
        dsmcStandardFields standardFields_;
        
        dsmcFieldProperties fields_; 

        //- boundaries
        dsmcBoundaries boundaries_;


        // flux-face tracking 
        dsmcFaceTracker trackingInfo_;

        // References to the cloud sub-models

        //- Binary collision model
        autoPtr<BinaryCollisionModel> binaryCollisionModel_;
        
        // - Collision selection routine model
        autoPtr<collisionPartnerSelection> collisionPartnerSelectionModel_;

        dsmcReactions reactions_;

        boundaryMeasurements boundaryMeas_;
        
    // Private Member Functions


        label pickFromCandidateList
        (
            DynamicList<label>& candidatesInCell
        );

        // whenever you pick from the candidate list 
        // always update the candidate sublist
        void updateCandidateSubList
        (
            const label& candidate,
            DynamicList<label>& candidatesInSubCell
        );

        label pickFromCandidateSubList
        (
            DynamicList<label>& candidatesInCell,
            DynamicList<label>& candidatesInSubCell
        );

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- Calculate collisions between molecules
        void collisions();

        //- Disallow default bitwise copy construct
        dsmcCloud(const dsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcCloud&);


public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //  file
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& dsmcInitialiseDict,
            const bool& clearFields
        );

/*
        //- Construct for deleteParcels
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            label dummy,
            bool readFields = true
        );

        //- Construct for mapping
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const vector& translationVector,
            const vector& scaleVector,
            bool readFields = true
        );*/

    //- Destructor
//     virtual ~dsmcCloud();


    //- Type of parcel the cloud was instantiated for
//     typedef dsmcParcel parcelType;


    // Member Functions

        // Access

    // References to the mesh and databases

        //- Return the cloud type
        inline const word& cloudName() const;

        //- Return refernce to the mesh
        inline const fvMesh& mesh() const;


    // References to the dsmc specific data

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;

        //- Return the idList
        inline const List<word>& typeIdList() const;

        //- Return the number of real particles represented by one
        //  parcel
        inline scalar nParticle() const;

        //- Return the cell occupancy addressing
        inline const List<DynamicList<dsmcParcel*> >&
        cellOccupancy() const;

        //- Return the sigmaTcRMax field.  non-const access to allow
        // updating.
        inline volScalarField& sigmaTcRMax();

        //- Return the collision selection remainder field.  non-const
        // access to allow updating.
        inline scalarField& collisionSelectionRemainder();

        //- Return all of the constant properties
        inline const List<dsmcParcel::constantProperties>&
        constProps() const;

        //- Return the constant properties of the given typeId
        inline const dsmcParcel::constantProperties&
        constProps(label typeId) const;

        //- Return refernce to the random object
        inline Random& rndGen();

    // Kinetic theory helper functions

        //- Generate a random velocity sampled from the Maxwellian speed
        // distribution
        vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

        //- Generate a random rotational energy, sampled from the
        // equilibrium distribution (Bird eqn 11.22 and 11.23 and
        // adapting code from DSMC3.FOR)
        scalar equipartitionRotationalEnergy
        (
            scalar temperature,
            scalar rotationalDegreesOfFreedom
        );
        
        //- Generate a random vibrational energy
        scalar equipartitionVibrationalEnergy
        (
            scalar temperature,
            scalar vibrationalDegreesOfFreedom,
            label typeId
        );

        scalar equipartitionVibrationalEnergy2
        (
            scalar temperature,
            scalar vibrationalDegreesOfFreedom,
            label typeId
        );

        vector chapmanEnskogVelocity
        (
            scalar temperature,
            scalar mass,
            vector q,
            tensor tau
        );

        vector chapmanEnskogVelocityMiu
        (
            scalar temperature,
            scalar mass,
            scalar B,
            vector q,
            tensor tau
        );

        void generalisedChapmanEnskog
        (
            label typeID,
            scalar translationalTemperature,
            scalar rotationalTemperature,
            scalar vibrationalTemperature,
            scalar mass,
            vector D,
            vector qTra,
            vector qRot,
            vector qVib,
            tensor tau,
            scalar& ERot,
            scalar& EVib,
            vector& U
        );

        void generalisedChapmanEnskog2
        (
            label typeID,
            scalar translationalTemperature,
            scalar rotationalTemperature,
            scalar vibrationalTemperature,
            scalar mass,
            vector D,
            vector qTra,
            vector qRot,
            vector qVib,
            tensor tau,
            scalar& ERot,
            scalar& EVib,
            vector& U
        );

        // From the Maxwellian distribution:
        //- Average particle speed
        inline scalar maxwellianAverageSpeed
        (
            scalar temperature,
            scalar mass
        ) const;

        inline scalarField maxwellianAverageSpeed
        (
            scalarField temperature,
            scalar mass
        ) const;

        //- RMS particle speed
        inline scalar maxwellianRMSSpeed
        (
            scalar temperature,
            scalar mass
        ) const;

        inline scalarField maxwellianRMSSpeed
        (
            scalarField temperature,
            scalar mass
        ) const;

        //- Most probable speed
        inline scalar maxwellianMostProbableSpeed
        (
            scalar temperature,
            scalar mass
        ) const;

        inline scalarField maxwellianMostProbableSpeed
        (
            scalarField temperature,
            scalar mass
        ) const;


    // Sub-models

        //- Return reference to binary elastic collision model
        inline const BinaryCollisionModel&  binaryCollision() const;

        //- Return non-const reference to binary elastic collision model
        inline BinaryCollisionModel& binaryCollision();

        inline const dsmcBoundaries& boundaries() const;  //****
        inline dsmcBoundaries& boundaries();//****

        inline const dsmcFieldProperties& fields() const;  //****
        inline dsmcFieldProperties& fields();//****

        inline const dsmcStandardFields& stdFields() const;  //****
        inline dsmcStandardFields& stdFields();//****

        inline const dsmcControllers& controllers() const;  //****
        inline dsmcControllers& controllers();//****

        inline const dsmcFaceTracker& tracker() const;  //****
        inline dsmcFaceTracker& tracker();//****
        
        inline const dsmcReactions& reactions() const;  //****
        inline dsmcReactions& reactions();//****

        
        inline const boundaryMeasurements& boundaryFluxMeasurements() const;
        inline boundaryMeasurements& boundaryFluxMeasurements();
        
        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational energy in the system
            inline scalar rotationalEnergyOfSystem() const;
        
        //- Total vibratational energy in the system
            inline scalar vibrationalEnergyOfSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;

            void insertParcelInCellOccupancy(dsmcParcel* p);

            void removeParcelFromCellOccupancy
            (
                const label& cellMolId,
                const label& cell
            );

            void reBuildCellOccupancy();





        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar ERot,
                const scalar EVib,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const label typeId,
                const label newParcel,
                const label classification
            );
            
            //- redistribute vibrational energy after a reaction
            scalar energyRatio
            (
                scalar ChiA,
                scalar ChiB
            );
            
            //- redistribute vibrational energy after a reaction
            scalar PSIm
            (
                scalar DOFm,
                scalar DOFtot
            );

            //- Evolve the cloud (move, collide)
            void evolve();
            void evolveWrite();

            void resetHybrid
            (
                volScalarField& TtrInitial,
                volVectorField& UInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& qInitial,
                PtrList<volTensorField>& tauInitial,
                dimensionedScalar& B,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybrid2
            (
                volScalarField& TtrInitial,
                volVectorField& UInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& qInitial,
                PtrList<volTensorField>& tauInitial,
                dimensionedScalar& B,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybridMax
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& numberDensitiesField,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybridTra
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& qtInitial,
                PtrList<volTensorField>& tauInitial,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybridTraRotVib
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& TrInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& DInitial,
                PtrList<volVectorField>& qtInitial,
                PtrList<volVectorField>& qrInitial,
                PtrList<volVectorField>& qvInitial,
                PtrList<volTensorField>& tauInitial,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybridTraRotVib2
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& TrInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& DInitial,
                PtrList<volVectorField>& qtInitial,
                PtrList<volVectorField>& qrInitial,
                PtrList<volVectorField>& qvInitial,
                PtrList<volTensorField>& tauInitial,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void shockReset();

            //- Clear the Cloud
            inline void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "dsmcCloud.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
