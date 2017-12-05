/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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


Description

   

SourceFiles
    polyFadeIII.H
    polyFadeII.C
    polyFadeIIIO.C

\*---------------------------------------------------------------------------*/

#ifndef polyFadeII_H
#define polyFadeII_H

#include "polyMoleculeCloud.H"
#include "IFstream.H"
#include "OFstream.H"
#include "cachedRandomMD.H"
#include "collisionDistribution.H"
#include "distributePoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyFadeII Declaration
\*---------------------------------------------------------------------------*/

class polyFadeII
{
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        cachedRandomMD rndGen_;

        label molId_;

        scalar deltaT_;
        
        scalar tauT_;

        label n_;
        
        label maxMolTries_;

        label nInserted_;
        label nDeleted_;

        bool boundedBox_;
        
        bool biggestRadius_;

        // insert
        bool closestPositionBoundBox_;
        bool closestPositionCell_;        
        bool exactPosition_;
        
        //delete
        bool deleteClosestPositionBoundBox_;
        bool deleteClosestPosition_;
        
        bool empty_;
        
        // stores information on all processors
        
        DynamicList<label> trackingNumbersIns_;
        DynamicList<label> trackingNumbersDel_;

        DynamicList<scalar> timeIns_;
        DynamicList<scalar> timeDel_;
        
        // insertion velocity
        scalar temperature_;
        vector velocity_;
        word velocityOption_;

        // outer list size of idList, inner list positions 
        // of existing molecules that are overlapping
        List< DynamicList<vector> > overlappingMolsPositions_;
        List< DynamicList<label> > overlappingMolsTNs_;        
        
    // Private Member Functions

        // insert

        polyMolecule* insertMoleculeInCloud
        (
            polyMoleculeCloud& molCloud,
            const vector& position,
            const label& cell,
            const tensor& Q,
            const vector& v,
            const vector& a,
            const vector& pi,
            const vector& tau,
            const vector& specialPosition,
            const label molId,
            const label special
        );
        
        void setVelocity
        (
            polyMoleculeCloud& molCloud,
            vector& velocity
        );
        
        vector getMaxwellianVelocity
        (
            polyMoleculeCloud& molCloud
        );        
        
        label insertMoleculeInBoundBox
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            DynamicList<polyMolecule*>& molsInBB,
            DynamicList<label>& failedTNs
        );
        
        label insertMoleculeAtClosestPositionBoundBox
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            DynamicList<polyMolecule*>& molsInBB,
            DynamicList<label>& failedTNs,
            const vector& r
        );
        
        label insertMoleculeAtClosestPositionCell
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            DynamicList<label>& failedTNs,
            const vector& r
        );        
        
        label insertMoleculeAtExactPosition
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            const vector& r
        );        
        
        // delete
        
        void deleteMoleculeFromCloud
        (
            polyMoleculeCloud& molCloud,
            polyMolecule& mol
        );
        
        label deleteMoleculeFromClosestPositionBoundBox
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            const vector& r
        );
        
        label deleteMoleculeFromClosestPosition
        (
            polyMoleculeCloud& molCloud,
            const boundedBox& bb,
            const vector& r
        );
        
        // updating
        
        void insertInLists
        (
            DynamicList<label>& tNsIns,
            DynamicList<label>& tNsDel,
            DynamicList<scalar>& timeIns,
            DynamicList<scalar>& timeDel 
        );

        
        void updateLists
        (
            DynamicList<label>& tNsIns,
            DynamicList<label>& tNsDel
        );
        
        // overlapping
        void deleteMolecules
        (
            polyMoleculeCloud& molCloud,
            const DynamicList<label>& tNs
        );       
        
        void readProperties();

        bool chooseRandomPoint
        (
            const polyMesh& mesh,
            const boundedBox& bb,
            const vector& rI, 
            const scalar& R,
            label& cellI,
            vector& rStart
        );
    
        void pickExistingMoleculeBiggestRadius
        (
            const DynamicList<polyMolecule*>& molsInBB,
            const DynamicList<label>& failedTNs,
            label& tN,
            vector& rI,
            scalar& R
        );
        
        void pickExistingMoleculeClosestPositionCell
        (
            polyMoleculeCloud& molCloud,
            const DynamicList<label>& failedTNs,
            label& tN,
            vector& rI,
            const vector& r,
            scalar& R
        );
        
        void pickExistingMoleculeClosestPositionBoundBox
        (
            const DynamicList<polyMolecule*>& molsInBB,
            const DynamicList<label>& failedTNs,
            label& tN,
            vector& rI,
            const vector& r,
            scalar& R     
        );          


public:

    // Static data members

    // Constructors

        //- Construct null
        polyFadeII();

        //- Construct from cachedRandomMD
		polyFadeII
		(
			const cachedRandomMD& rndMD
		);

        //- Construct from dict
        polyFadeII
        (
            const dictionary& dict,
			const cachedRandomMD& rndMD
        );


    // Selectors

    // Destructor

        ~polyFadeII();

    // Member Functions

        void createInitialConfiguration
        (
            const dictionary& dict,        
            const label& molId,
            const scalar& deltaT,
            const word& insertOption,
            const word& deleteOption
        );
        
        void setDeleteOption(const word& option);
        
        void setInsertOption(const word& option);        

        void controlMolecules
        (
            polyMoleculeCloud& molCloud,
            const label& insertOrDelete,
            const boundedBox& bb,
            const List<vector>& positions 
        );


        void updateProperties
        (
            const dictionary& newDict
        );

        void checkFractions(polyMoleculeCloud& molCloud);
        
        // Access

            inline const label& nInserted() const;
            inline const label& nDeleted() const;
            inline const DynamicList<label>& insertingList() const;
            inline const DynamicList<label>& deletingList() const;
            inline const List< DynamicList<vector> >& overlappingMols() const;
            inline bool isNotFading(const label& tN);
            inline const scalar& temperature() const;
            inline scalar& temperature();
            inline const vector& velocity() const;
            inline vector& velocity();
            
        // Check

        // Edit
            
        // Write

            void write
            (
                const Time& time,
                const fileName& fixedPathName,
                const fileName& timePath
            );
            
            void output(Time& time);

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyFadeIII.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
