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
    Foam::polyMoleculeCloud

Description

SourceFiles
    polyMoleculeCloudI.H
    polyMoleculeCloud.C


\*---------------------------------------------------------------------------*/

#ifndef polyMoleculeCloud_H
#define polyMoleculeCloud_H

#include "Cloud.H"
#include "polyMolecule.H"
#include "IOdictionary.H"
#include "labelVector.H"
#include "cachedRandomMD.H"
#include "fileName.H"
#include "reducedUnits.H"
#include "polyFieldProperties.H"
#include "polyBoundaries.H"
#include "polyControllers.H"
#include "polyFaceTracker.H"
#include "trackingNumber.H"
#include "cyclicBoundaries.H"
#include "cellInteractions.H"
#include "clockTimer.H"
#include "potentials.H" 
#include "constantMoleculeProperties.H"
#include "polyIntegrators.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyMoleculeCloud Declaration
\*---------------------------------------------------------------------------*/

class polyMoleculeCloud
:
    public Cloud<polyMolecule>
{

private:

    // Private data

        const polyMesh& mesh_;
        
        const reducedUnits& redUnits_;

        const constantMoleculeProperties& cP_;
        
        cachedRandomMD& rndGen_;

        polyIntegrators int_;
        
        potentials p_;
        
        List<DynamicList<polyMolecule*> > cellOccupancy_;
        
        //- measurement fields
        polyFieldProperties fields_; 
        
        
        //- boundaries
        polyBoundaries boundaries_;
        
        //- controllers
        polyControllers controllers_;
        
        //- polyMolecule tracking
        polyFaceTracker trackingInfo_;
        
        trackingNumber moleculeTracking_;
        
        cyclicBoundaries cyclics_;
        
        cellInteractions<polyMolecule> iL_;
        
        //- interaction particle list
        List<DynamicList<polyMolecule*> > ipl_;
        
        clockTimer clock_;


    // Private Member Functions

        void setSiteSizesAndPositions();
        void calculatePairForces();
        void checkForOverlaps();
        void checkMoleculesInMesh();
        void setMoleculesForDirectInteraction();
        void computePairForces();
        void accelerationUpdate();
        void velocityUpdate(const scalar& trackTime);


        //- Disallow default bitwise copy construct
        polyMoleculeCloud(const polyMoleculeCloud&);
        //- Disallow default bitwise assignment
        void operator=(const polyMoleculeCloud&);

public:

    // Constructors

        //- Construct for mdFoam
        polyMoleculeCloud
        (
            Time& t,
            const polyMesh& mesh,
//             const potentials& p,
            const reducedUnits& rU,
            const constantMoleculeProperties& cP,
            cachedRandomMD& rndGen
        );

        //- General constructor
        polyMoleculeCloud
        (
            Time& t,
            const polyMesh& mesh,
//             const potentials& p,
            const reducedUnits& rU,
            const constantMoleculeProperties& cP,
            cachedRandomMD& rndGen,
            const word& option,
            const bool& clearFields
        );
        
    // Static constructors

        static autoPtr<polyMoleculeCloud> New
        (
            Time& t,
            const polyMesh& mesh,
//             const potentials& p,
            const reducedUnits& rU,
            const constantMoleculeProperties& cP,
            cachedRandomMD& rndGen
        );

        static autoPtr<polyMoleculeCloud> New
        (
            Time& t,
            const polyMesh& mesh,
//             const potentials& p,
            const reducedUnits& rU,
            const constantMoleculeProperties& cP,
            cachedRandomMD& rndGen,
            const word& option,
            const bool& clearFields
        );


    // Member Functions

        void createMolecule
        (
            const vector& position,
            const label cell,
            const label tetFace,
            const label tetPt,         
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
            const label special,
            const label id,
            const scalar& fraction,
            const label trackingNumber
        );
        
        void updateTrackingNumbersAfterRead();
        
        //- Evolve the molecules (move, calculate forces, etc)
        
        void evolve();
        void evolveBeforeForces();
        void evolveAfterForces();
        
        void controlBeforeVelocity();
        void updateVelocity();
        void controlBeforeMove();
        void move();
        void move(const scalar& trackTime);        
        void updateAfterMove(const scalar& trackTime);        
        void controlAfterMove();
        void controlBeforeForces();
        void clearLagrangianFields();
        void calculateForce();
        void updateAcceleration();
        void controlAfterForces();
        void controlAfterVelocity();
        void postTimeStep();

        
        void prepareInteractions();
        void buildCellOccupancy();
        void rebuildCellOccupancy();

        inline void evaluatePair
        (
            polyMolecule* molI,
            polyMolecule* molJ
        );

        inline bool evaluatePotentialLimit
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            const scalar& potentialEnergyLimit
        ) ;

        void insertMolInCellOccupancy
        (
            polyMolecule* mol
        );

        void removeMolFromCellOccupancy
        (
            const label& cellMolId,
            const label& cell
        );

        void removeMolFromCellOccupancy
        (
            polyMolecule* molI
        );

        void updateNeighbouringRadii
        (
            polyMolecule* molI
        );

        void updateRadii();
        void readNewField();
        void setIPL();
        label getTrackingNumber();
        void resetTrackingNumbers();
        void writeReferredCloud();
        
        // Access

            label nSites() const;
            inline const polyMesh& mesh() const;
            inline const potentials& pot() const;
            inline potentials& pot();            
            inline const constantMoleculeProperties& cP() const;
//             inline const potential& pot() const;
            inline const List<DynamicList<polyMolecule*> >& cellOccupancy() const;
//             inline const List<polyMolecule::constantProperties> constProps() const;
//             inline const polyMolecule::constantProperties& constProps(label id) const;
            inline cachedRandomMD& rndGen();
            inline const reducedUnits& redUnits() const;
            inline const polyFieldProperties& fields() const;  
            inline polyFieldProperties& fields();
            inline const cyclicBoundaries& cyclics() const;
            inline const polyBoundaries& boundaries() const;
            inline polyBoundaries& boundaries();
            inline const polyControllers& controllersInfo() const; 
            inline polyControllers& controllersInfo();
            inline const polyFaceTracker& tracker() const;
            inline polyFaceTracker& tracker();
            inline trackingNumber& moleculeTracking();
            inline cellInteractions<polyMolecule>& il();
            inline clockTimer& clock();
            
    // Member Operators

        //- Write polyMolecule sites in XYZ format
        void writeXYZ(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMoleculeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
