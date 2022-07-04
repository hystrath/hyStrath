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
    Foam::pdCloud

Description
    Templated base class for pd cloud

SourceFiles
    pdCloudI.H
    pdCloud.C

\*---------------------------------------------------------------------------*/

#ifndef pdCloud_H
#define pdCloud_H

#include "Cloud.H"

#include "pdParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"

#include "volFields.H"
#include "scalarIOField.H"
#include "BinaryCollisionModel.H"

#include "pdBoundaries.H"
#include "pdAllConfigurations.H"
#include "pdFieldProperties.H"
#include "pdFaceTracker.H"
#include "pdControllers.H"
#include "pdStandardFields.H"
/***********************************************************************/
#include "pdEmFields.H"
#include "ElectronModel.H"
//#include "nAveragingMethod.H"
/***********************************************************************/

#include "pdReactions.H"
#include "collisionPartnerSelection.H"
// #include "reactInfoList.H"
#include "boundaryMeasurements.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pdCloud Declaration
\*---------------------------------------------------------------------------*/


class pdCloud
:
    public Cloud<pdParcel>
{

private:

    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word& cloudName_;

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
        List<DynamicList<pdParcel*> > cellOccupancy_;

        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<pdParcel::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- controllers
        pdControllers controllers_;

        //- measurement fields
        pdStandardFields standardFields_;

        /**************************************************************/
        //- electromagnetic fields
        pdEmFields emFields_;

        //- Electron fluid model
        autoPtr<ElectronModel> electronModel_;
        /**************************************************************/

        pdFieldProperties fields_;

        //- boundaries
        pdBoundaries boundaries_;


        // flux-face tracking
        pdFaceTracker trackingInfo_;

        // References to the cloud sub-models

        //- Binary collision model
        autoPtr<BinaryCollisionModel> binaryCollisionModel_;

        // - Collision selection routine model
        autoPtr<collisionPartnerSelection> collisionPartnerSelectionModel_;

        pdReactions reactions_;

        boundaryMeasurements boundaryMeas_;

        const dictionary measureDict_;

        bool    stdMeasure_;
        label   stdFrequency_;
        label   stdCount_;

        label   infoFrequency_;


    // Private Member Functions


        label pdkFromCandidateList
        (
            DynamicList<label>& candidatesInCell
        );

        // whenever you pdk from the candidate list
        // always update the candidate sublist
        void updateCandidateSubList
        (
            const label& candidate,
            DynamicList<label>& candidatesInSubCell
        );

        label pdkFromCandidateSubList
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
        pdCloud(const pdCloud&);

        //- Disallow default bitwise assignment
        void operator=(const pdCloud&);


public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //  file
        pdCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        pdCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& pdInitialiseDict,
            const bool& clearFields
        );

/*
        //- Construct for deleteParcels
        pdCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            label dummy,
            bool readFields = true
        );

        //- Construct for mapping
        pdCloud
        (
            Time& t,
            const word& cloudName,
            const fvMesh& mesh,
            const vector& translationVector,
            const vector& scaleVector,
            bool readFields = true
        );*/

    //- Destructor
     //~pdCloud();

    // Member Functions

        // Access

    // References to the mesh and databases

        //- Return the cloud type
        inline const word& cloudName() const;

        //- Return refernce to the mesh
        inline const fvMesh& mesh() const;


    // References to the pd specific data

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;

        //- Return the idList
        inline const List<word>& typeIdList() const;

        //- Return the number of real particles represented by one
        //  parcel
        inline scalar nParticle() const;

        //- return info output frequency
        inline label infoFrequency() const;

        //- Return the cell occupancy addressing
        inline const List<DynamicList<pdParcel*> >&
        cellOccupancy() const;

        //- Return the sigmaTcRMax field.  non-const access to allow
        // updating.
        inline volScalarField& sigmaTcRMax();

        //- Return the collision selection remainder field.  non-const
        // access to allow updating.
        inline scalarField& collisionSelectionRemainder();

        //- Return all of the constant properties
        inline const List<pdParcel::constantProperties>&
        constProps() const;

        //- Return the constant properties of the given typeId
        inline const pdParcel::constantProperties&
        constProps(label typeId) const;

        //- Return refernce to the random object
        inline Random& rndGen();

    // Random implementation corrections
        //Generates a random integer between valMin and valMax
        label randomLabel
        (
            const label valOne,
            const label valTwo
        );

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
        // adapting code from PD3.FOR)
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

        inline const pdBoundaries& boundaries() const;  //****
        inline pdBoundaries& boundaries();//****

        inline const pdFieldProperties& fields() const;  //****
        inline pdFieldProperties& fields();//****

        inline const pdStandardFields& stdFields() const;  //****
        inline pdStandardFields& stdFields();//****

        inline const pdEmFields& emFields() const;  //****
        inline pdEmFields& emFields();//****

        //- Return reference to binary elastic collision model
        inline const ElectronModel&  electronModel() const;
         //- Return non-const reference to electron fluid model
        inline ElectronModel& electronModel();
        /*
        //- Return reference to binary elastic collision model
        inline const AveragingMethod&  averagingMethod() const;
         //- Return non-const reference to electron fluid model
        inline AveragingMethod& averagingMethod();
        */

        inline const pdControllers& controllers() const;  //****
        inline pdControllers& controllers();//****

        inline const pdFaceTracker& tracker() const;  //****
        inline pdFaceTracker& tracker();//****

        inline const pdReactions& reactions() const;  //****
        inline pdReactions& reactions();//****


        inline const boundaryMeasurements& boundaryFluxMeasurements() const;
        inline boundaryMeasurements& boundaryFluxMeasurements();

        // Check

            //- Total mass injected
            //inline scalar massInjected() const;

            //- Total mass in system
            //inline scalar massInSystem() const;

            //- Total mass in system
            //inline scalar chargeInSystem() const;

            //- Total linear momentum of the system
            //inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            //inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational energy in the system
            //inline scalar rotationalEnergyOfSystem() const;

        //- Total vibratational energy in the system
            //inline scalar vibrationalEnergyOfSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;

            void insertParcelInCellOccupancy(pdParcel* p);

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
                const vector& A,
                const scalar EPot,
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

            //- setup PIC method
            void setupPIC();

            //- Clear the Cloud
            inline void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "pdCloud.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
