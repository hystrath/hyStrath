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
    Foam::pdParcel

Description
    PD parcel class

SourceFiles
    pdParcelI.H
    pdParcel.C
    pdParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef pdParcel_H
#define pdParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

// #include "PdCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// template<class pdParcel>
// class pdParcel;

class pdCloud;

// Forward declaration of friend functions

// template<class pdParcel>
// Ostream& operator<<
// (
//     Ostream&,
//     const pdParcel<pdParcel>&
// );

/*---------------------------------------------------------------------------*\
                         Class pdParcel Declaration
\*---------------------------------------------------------------------------*/

// template<class pdParcel>
class pdParcel
:
//     public Particle<pdParcel>
      public particle
{
public:

    //- Class to hold PD particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle mass [kg]                  (constant)
            scalar mass_;

            //- Particle charge units []            (constant)
            scalar Ze_;

            //- Particle hard sphere diameter [m]   (constant)
            scalar d_;

            //- Rotational degrees of freedom
            scalar rotationalDegreesOfFreedom_;

            //- Vibrational degrees of freedom
            scalar vibrationalDegreesOfFreedom_;

            //- Viscosity index
            scalar omega_;

            //- VSS scattering parameter
            scalar alpha_;

            //- Characteristic vibrational temperature
            scalar thetaV_;

            //- characteristic dissociation quantum level
            label charDissQuantumLevel_;

            //- Characteristic dissociation temperature
            scalar thetaD_;

            //- Reference vibrational relaxation number
            scalar Zref_; //TODO: Change label to avoid confusion with Ze_

            //- Reference temperature for reference vibrational relaxation number
            scalar TrefZv_;

            label newParcel_;

            label classification_;


        public:

        // Constructors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);

            //inline ~constantProperties();

        // Member functions

            //- Return const access to the particle density
            inline scalar mass() const;

            //- Return const access to the particle charge units
            inline scalar Ze() const;

            //- Return const access to the minimum particle mass
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the rotationalDegreesOfFreedom
            inline scalar rotationalDegreesOfFreedom() const;

            //- Return the vibrationallDegreesOfFreedom
            inline scalar vibrationalDegreesOfFreedom() const;

            inline label charDissQuantumLevel() const;

            //- Return the viscosity index
            inline scalar omega() const;

            //- Return the VSS scattering parameter
            inline scalar alpha() const;

            //- Return the characteristic vibrational temperature
            inline scalar thetaV() const;

            //- Return the characteristic dissociation temperature
            inline scalar thetaD() const;

            //- Return the characteristic vibrational relaxation number
            inline scalar Zref() const;

            //- Return the reference temperature for the characteristic vibrational relaxation number
            inline scalar TrefZv() const;

    };

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<pdCloud>
    {
        //Private data

        // label specifying which part of the integration algorithm is taking place
        label part_;

        public:

        // Constructors
            //- Construct from components
            inline trackingData(pdCloud& cloud, label part);

            /*
            trackingData(pdCloud& cloud, label part)
            :
                particle::TrackingData<pdCloud>(cloud),
                part_(part)
            {}
            */
        // Member functions

            //- Return leapfrog integration tag
            inline label part() const;

    };

    private:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Acceleration of Parcel [m/s^2]
            vector A_;

            //- Potential energy of the Parcel [V]
            scalar EPot_;

            //- Rotational energy of the Parcel [J]
            scalar ERot_;

            //- Vibrational energy of the Parcel [J]
            scalar EVib_;

            //- Parcel type id
            label typeId_;

            //- newly inserted parcel
            label newParcel_;

            // see W.F. Santos, "Physical and computational aspects
            // of shock waves over power-law leading edges",
            // Physics of Fluids, 20, 016101, 2008.
            label classification_;


public:

    friend class Cloud<pdParcel>;


    // Constructors

        //- Construct from components
        inline pdParcel
        (
            const polyMesh& mesh,
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

        pdParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new pdParcel(*this));
        }


        //- Factory class to read-construct particles used for parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<pdParcel> operator()(Istream& is) const
            {
                return autoPtr<pdParcel>(new pdParcel(mesh_, is, true));
            }
        };


    // Destructor

        //~pdParcel();

    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to acceleration
            inline const vector& A() const;

            //- Return const access to potential energy
            inline scalar EPot() const;

            //- Return const access to rotational energy
            inline scalar ERot() const;

            //- Return const access to vibrational energy
            inline scalar EVib() const;

            //- Return const access to new parcel
            inline label newParcel() const;

            //- Return const access to parcel classification
            // - see Santos W.F., Physical and computational
            //- aspects of shock waves over power-law leading edges,
            //- Physics of Fluids, 20, 016101, 2008
            inline label classification() const;

        // Edit

            //- Return non-const access to typeId
            inline label& typeId();

            //- Return access to velocity
            inline vector& U();

            //- Return access to acceleration
            inline vector& A();

            //- Return access to potential energy
            inline scalar& EPot();

            //- Return access to rotational energy
            inline scalar& ERot();

            //- Return access to vibrational energy
            inline scalar& EVib();

            //- Return access to new parcel
            inline label& newParcel();

            //- Return access to parcel classification
            // - see Santos W.F., Physical and computational
            //- aspects of shock waves over power-law leading edges,
            //- Physics of Fluids, 20, 016101, 2008
            inline label& classification();


        // Main calculation loop

            // Tracking

                //- Move the parcel
                //template<class trackData>
            bool move(trackingData& td, const scalar trackTime);

        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            //template<class trackData>
            // Patch interactions

            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);


        // I-O

            static void readFields(Cloud<pdParcel>& c);

            static void writeFields(const Cloud<pdParcel>& c);


    // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream&,
            const pdParcel&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pdParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
