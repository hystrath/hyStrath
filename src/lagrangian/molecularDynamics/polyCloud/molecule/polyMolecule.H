/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::polyMolecule

Description
    Foam::polyMolecule

SourceFiles
    polyMoleculeI.H
    polyMolecule.C
    polyMoleculeIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyMolecule_H
#define polyMolecule_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
// #include "diagTensor.H"
// #include "reducedUnits.H"
#include "constantMoleculeProperties.H"

// #include "constPropSite.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMoleculeCloud;

/*---------------------------------------------------------------------------*\
                           Class polyMolecule Declaration
\*---------------------------------------------------------------------------*/

class polyMolecule
:
    public particle
{

public:

    // Values of special that are less than zero are for built-in functionality.
    // Values greater than zero are user specifiable/expandable (i.e. test
    // special_ >= SPECIAL_USER)

    enum specialTypes
    {
        SPECIAL_TETHERED = -1,
        SPECIAL_FROZEN   = -2,
        NOT_SPECIAL      = 0,
        SPECIAL_REFERRED = 1,
        SPECIAL_USER     = 2
    };

    /*
    //- Class to hold molecule constant properties
    class constantProperties
    {

        // Private data

            label nSites_;

            List<constPropSite> sites_;

            List<label> electrostaticSites_;

            List<label> pairPotSites_;

            diagTensor momentOfInertia_;

            scalar mass_;

        // Private Member Functions

            void checkSiteListSizes() const;

            void setInteracionSiteBools
            (
                const List<word>& siteIds,
                const List<word>& pairPotSiteIds
            );

            bool linearMoleculeTest() const;


    public:

        inline constantProperties();

        //- Construct from dictionary
        inline constantProperties
        (
            const dictionary& dict,
            const reducedUnits& rU,
            const labelList& siteIds
        );

        // Member functions

        inline const List<constPropSite>& sites() const;

        inline const List<label>& pairPotSites() const;

        inline const List<label>& electrostaticSites() const;

        inline const diagTensor& momentOfInertia() const;

        inline bool linearMolecule() const;

        inline bool pointMolecule() const;

        inline label degreesOfFreedom() const;

        inline scalar mass() const;

        inline label nSites() const;
    };*/


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<polyMoleculeCloud>
    {
        // label specifying which part of the integration algorithm is taking
        label part_;


    public:

        // Constructors

            trackingData(polyMoleculeCloud& cloud, label part)
            :
                particle::TrackingData<polyMoleculeCloud>(cloud),
                part_(part)
            {}

        // Member functions

            inline label part() const
            {
                return part_;
            }
    };    
    

private:

    // Private data

        //- Be careful with the ordering of data.
        //  It has an impact on binary transfer:
        //    -# Put the largest data members 1st
        //    -# Pair up labels,
        //    -# Don't go scalar-label, scalar-label, because in 64bit mode,
        //       the labels will be padded by 4bytes.

        tensor Q_;

        //- linear velocity
        vector v_;

        //- linear acceleration
        vector a_;

        //- angular momentum
        vector pi_;

        //- torque
        vector tau_;

        vector specialPosition_;

        scalar potentialEnergy_;

        // - r_ij f_ij, stress dyad
        tensor rf_;

        label special_;

        label id_;

        scalar R_;

        // - fraction of force potential contribution
        scalar frac_;    

        //- each molecule will have a unique label
        label trackingNumber_;

        List<vector> siteForces_;

        List<vector> sitePositions_;


    // Private Member Functions

        tensor rotationTensorX(scalar deltaT) const;

        tensor rotationTensorY(scalar deltaT) const;

        tensor rotationTensorZ(scalar deltaT) const;


public:

    friend class Cloud<polyMolecule>;

    // Constructors

        //- Construct from components
        inline polyMolecule
        (
            const polyMesh& mesh,            
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
//             const constantProperties& constProps,
            const constantMoleculeProperties& cP,
            const label special,
            const label id,
            const scalar& frac,
            const label trackingNumber
        );

        //- Construct from Istream
        polyMolecule
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new polyMolecule(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<polyMolecule> operator()(Istream& is) const
            {
                return autoPtr<polyMolecule>(new polyMolecule(mesh_, is, true));
            }
        };

    // Member Functions

        // Tracking

            bool move(trackingData&, const scalar& trackTime);

            void transformProperties(const tensor& T);

            void transformProperties(const vector& separation);

            void setSitePositions(const constantMoleculeProperties& cP);

            void setSiteSizes(label size);

            void updateHalfVelocity
            (
                const constantMoleculeProperties& cP,
                const scalar& trackTime
            );

            void updateAcceleration
            (
                const constantMoleculeProperties& cP
            );

            void updateAfterMove
            (
                const constantMoleculeProperties& cP,
                const scalar& trackTime
            );

            void setAsReferred();
            
        // Access

            inline const tensor& Q() const;
            inline tensor& Q();

            inline const vector& v() const;
            inline vector& v();

            inline const vector& a() const;
            inline vector& a();

            inline const vector& pi() const;
            inline vector& pi();

            inline const vector& tau() const;
            inline vector& tau();

            inline const List<vector>& siteForces() const;
            inline List<vector>& siteForces();

            inline const List<vector>& sitePositions() const;
            inline List<vector>& sitePositions();

            inline const vector& specialPosition() const;
            inline vector& specialPosition();

            inline scalar potentialEnergy() const;
            inline scalar& potentialEnergy();

            inline const tensor& rf() const;
            inline tensor& rf();

            inline label special() const;
            inline label& special();
            
            inline bool tethered() const;

            inline label id() const;
            inline label& id();

            // - return shortest radius
            inline const scalar& R() const;
            inline scalar& R();

            inline const scalar& fraction() const;
            inline scalar& fraction();

            // - return tracking number
            inline const label& trackingNumber() const;
            inline label& trackingNumber();

            inline bool frozen();
            
            inline bool referred();
            
    // Member Operators


        //- Overridable function to handle the particle hitting a patch
        //  Executed before other patch-hitting functions
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

        static void readFields(Cloud<polyMolecule>& mC);

        static void writeFields(const Cloud<polyMolecule>& mC);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const polyMolecule&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMoleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
