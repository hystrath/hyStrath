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
#include "dynamicFvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "BinaryCollisionModel.H"

#include "dsmcBoundaries.H"
#include "dsmcAllConfigurations.H"
#include "dsmcFieldProperties.H"
#include "dsmcFaceTracker.H"
#include "dsmcControllers.H"
#include "dsmcDynamicLoadBalancing.H"

#include "dsmcReactions.H"
#include "collisionPartnerSelection.H"
#include "boundaryMeasurements.H"
#include "cellMeasurements.H"

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
        const word& cloudName_;

        //- References to the mesh and time databases
        const dynamicFvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;
        
        //- Dictionary of particle properties
        IOdictionary controlDict_;

        //- A list of unique instances of molecule types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Employ the variable time step method?
        bool variableTimeStepMethod_;
        
        //- Number of real atoms/molecules represented by a parcel
        //scalar nParticle_;
        volScalarField nParticle_;
        
        //- Variable time step method
        volScalarField variableTimeStep_;
        
        //- Essential info to store for the Variable time step method
        //  when mesh has changed: [0]: ref volume, [1]: ref number density
        //  [2]: ref deltaT
        FixedList<scalar, 3> refQuantitiesVariableTimeStep_;
        
        //- Bool to identify is this is an axisymmetric simulation
        //- Movements dealt with differently if true
        bool axisymmetric_;
        
        //- Radial extent of the domain, needed for axisymmetric
        //- simulations
        scalar radialExtent_;
        
        //- Maximum radial scaling factor, needed for axisymmetric
        //- simulations
        scalar maxRWF_;
        
        //- Write out the information to terminal every nTerminalOutputs
        label nTerminalOutputs_;

        //- A data structure holding which particles are in which cell
        //List<DynamicList<dsmcParcel*> > cellOccupancy_;
        DynamicList<DynamicList<dsmcParcel*> > cellOccupancy_; // EDITED VINCENT
        
        //- A data structure holding electron number density
        scalarList rhoNMeanElectron_;
        
        //- A data structure holding electron mass density
        scalarList rhoMMeanElectron_;
        
        //- A data structure holding all species mass density
        scalarList rhoMMean_;
        
        //- A data structure holding electron momentum
        List<vector> momentumMeanElectron_;
        
        //- A data structure holding all species momentum
        List<vector> momentumMean_;
        
        //- A data structure holding electron linearKE
        scalarList linearKEMeanElectron_;
        
        //- A data structure holding electron temperature
        scalarList electronTemperature_;
        
        //- A data structure holding cell velocity
        List<vector> cellVelocity_;

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
        
        //- load balancing
        dsmcDynamicLoadBalancing dynamicLoadBalancing_;

        //- measurement fields
        //dsmcStandardFields standardFields_;
        dsmcFieldProperties fields_; 

        //- boundaries
        dsmcBoundaries boundaries_;

        //- flux-face tracking 
        dsmcFaceTracker trackingInfo_;


    // References to the cloud sub-models

        //- Binary collision model
        autoPtr<BinaryCollisionModel> binaryCollisionModel_;
        
        // - Collision selection routine model
        autoPtr<collisionPartnerSelection> collisionPartnerSelectionModel_;

        dsmcReactions reactions_;

        boundaryMeasurements boundaryMeas_;
        
        //used to calculate mean collision separation
        cellMeasurements cellMeas_;
        
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

        //- When mesh has changed, changed the size of cellOccupancy_
        void buildCellOccupancyFromScratch();
        
        //- When mesh has changed, change the size of collisionSelectionRemainder_
        void buildCollisionSelectionRemainderFromScratch();
        
        //- When mesh has changed, update the boundaries
        void resetBoundaries();
        
        //- When mesh has changed, update the measurement tools
        void resetMeasurementTools();
        
        //- Remove electrons after adding their velocities to a DynamicList
        void removeElectrons();
        
        //- Add an electron where an ion is found, after the move function
        void addElectrons();

        //- Calculate collisions between molecules
        void collisions();
        
        //- Calculate number of equivalent particles field according to the
        //  variable time step method described in [Wu, 2004]
        void calculatenParticleField();
        
        //- Set the variable time step field according to the variable
        //  time-step method described in [Wu, 2004]
        //void setVariableTimeStepField(const label refCell);
        
        //- Calculate the variable time step field according to the 
        //  variable time-step method described in [Wu, 2004]
        void calculateTimeStepField();
        
        //- Update time step and nParticle fields
        //  when the mesh is changed
        void updateVariableTimeStep();

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
            const dynamicFvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const dynamicFvMesh& mesh,
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
     virtual ~dsmcCloud();


    //- Type of parcel the cloud was instantiated for
//     typedef dsmcParcel parcelType;


    // Member Functions

        // Access

    // References to the mesh and databases

        //- Return the cloud type
        inline const word& cloudName() const;

        //- Return refernce to the mesh
        inline const dynamicFvMesh& mesh() const;


    // References to the dsmc specific data

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;
        
        //- Return particle properties dictionary
        inline const IOdictionary& controlDict() const;

        //- Return the idList
        inline const List<word>& typeIdList() const;

        //- Return the number of real particles represented by one
        //  parcel
        //inline scalar nParticle() const;
        inline volScalarField& nParticle();
        
        //- Return the number of real particles represented by one
        //  parcel in a given cell
        inline scalar nParticle(const label celli);
        
        //- Return the const variable time step field
        inline const scalarField& deltaTValue();
        
        //- Return the time step value in a given cell
        inline scalar deltaTValue(const label celli);
        
        //- Return the radius of the domain for axisymmetric cases
        inline scalar radialExtent() const;
        
        //- Return the RWF for axisymmetric cases
        inline scalar maxRWF() const;
        
         //- Return the bool to identify if this is an axisymmetric simulation
        //- Measurements, etc dealt with differently if true
        inline bool axisymmetric() const;
        
        //- Return the number of timesteps to write out to terminal at
        inline label nTerminalOutputs() const;

        //- Return the cell occupancy addressing
        inline const DynamicList<DynamicList<dsmcParcel*> >& //EDITED VINCENT
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
        labelList equipartitionVibrationalEnergyLevel
        (
            scalar temperature,
            scalar vibrationalDegreesOfFreedom,
            label typeId
        );
        
        //- Generate a random electronic energy level
        label equipartitionElectronicLevel
        (
            scalar temperature,
            List<label> degeneracyList_,
            List<scalar> electronicEnergyList_,
            label typeId
        );
        
        // - Generate a LB post-collision rotational energy
        scalar postCollisionRotationalEnergy
        (
            scalar rotationalDof,
            scalar ChiB
        );
        
        // - Generate a LB post-collision vibrational energy
        label postCollisionVibrationalEnergyLevel
        (
            bool postReaction,
            label vibLevel,
            label iMax,
            scalar thetaV,
            scalar thetaD,
            scalar refTempZv,
            scalar omega,
            scalar Zref,
            scalar Ec
        );
        
        // - Generate a LB post-collision electronic energy level
        label postCollisionElectronicEnergyLevel
        (
            scalar Ec,
            label jMax,
            scalar omega,
            List<scalar> EElist,
            List<label> gList
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

        inline const dsmcControllers& controllers() const;  //****
        inline dsmcControllers& controllers();//****
        
        inline const dsmcDynamicLoadBalancing& dynamicLoadBalancing() const;  //****
        inline dsmcDynamicLoadBalancing& dynamicLoadBalancing();//****

        inline const dsmcFaceTracker& tracker() const;  //****
        inline dsmcFaceTracker& tracker();//****
        
        inline const dsmcReactions& reactions() const;  //****
        inline dsmcReactions& reactions();//****

        
        inline const boundaryMeasurements& boundaryFluxMeasurements() const;
        inline boundaryMeasurements& boundaryFluxMeasurements();
        
        inline const cellMeasurements& cellPropMeasurements() const;
        inline cellMeasurements& cellPropMeasurements();
        
        // Check

            //- Total mass in system
            inline scalarList infoMeasurements() const;

            //- Print cloud information
            void info() const;
            
             //- Check if load balance is required
            void loadBalanceCheck();
            
            //- Perform load balance
            void loadBalance(const int noRefinement = 0);

	// Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            void autoMap(const mapPolyMesh&);


            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;

            void insertParcelInCellOccupancy(dsmcParcel* p);

            void removeParcelFromCellOccupancy
            (
                const label& cellMolId,
                const label& cell
            );

            void reBuildCellOccupancy();
            
            void axisymmetricWeighting();


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar RWF,
                const scalar ERot,
                const label ELevel,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const label typeId,
                const label newParcel,
                const label classification,
                const labelList vibLevel
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
            
            // Retrun how often to write to terminal
            // (used for the solver to access this)
            label nTerminalOutputs();

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
