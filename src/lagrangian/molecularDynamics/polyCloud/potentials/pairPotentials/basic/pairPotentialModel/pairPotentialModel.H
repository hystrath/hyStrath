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
    pairPotentialModel

Description

SourceFiles
    pairPotentialModelI.H
    pairPotentialModel.C
    pairPotentialModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef pairPotentialModel_H
#define pairPotentialModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "writeTimeData.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"
#include "reducedUnits.H"
#include "polyMolecule.H"
#include "exclusionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
class energyScalingFunction;

/*---------------------------------------------------------------------------*\
                           Class pairPotentialModel Declaration
\*---------------------------------------------------------------------------*/

class pairPotentialModel
{

protected:

    // Protected data

        const polyMesh& mesh_;

        polyMoleculeCloud& molCloud_;
        
        const reducedUnits& rU_;
        
        dictionary  pairPotentialProperties_;
       
        // ids header of dictionary e.g. A-B 
        word name_;
        
        // split ids [0] = A; [1] = B
        List<word> idList_;
        
        scalar rCut_;
        scalar rCutSqr_;
        
        scalar rMin_;
        
        scalar dr_;        
        
        bool useTables_;
        
//         word energyScalingFunction_;


        scalar fMin_;
        scalar energyMin_;        
       
        List<scalar> forceLookup_;
        List<scalar> energyLookup_;
        
        mutable energyScalingFunction* esfPtr_;

        bool writeTables_;        

        bool exclusions_;        
        
        autoPtr<exclusionModel> exclusionModel_;
        
        
    // Private Member Functions
        
        void scaleEnergy
        (
            scalar& e,
            const scalar r
        ) const;        

public:

    //- Runtime type information
    TypeName("pairPotentialModel");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            pairPotentialModel,
            dictionary,
            (
                const polyMesh& mesh,
                polyMoleculeCloud& molCloud,
                const reducedUnits& redUnits,
                const word& name,
                const dictionary& dict
            ),
            (mesh, molCloud, redUnits, name, dict)
        );
   
    // Constructors

        //- Construct from components
        pairPotentialModel
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const reducedUnits& redUnits,
            const word& name,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<pairPotentialModel> New
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const reducedUnits& redUnits,
            const word& name,
            const dictionary& dict
        );



    // Destructor

        virtual ~pairPotentialModel();


    // Member Functions
        
        void setLookupTables();



        scalar energyLookUpFromTable (const scalar r) const;
        scalar forceLookUpFromTable (const scalar r) const;
        
        virtual scalar energy (const scalar r) const = 0;
        virtual scalar force (const scalar r) const = 0;
//         virtual scalar interaction (const scalar r, scalar& force, scalar& energy) const = 0;        
        
        
        List<Pair<scalar> > energyTable() const;
        List<Pair<scalar> > forceTable() const;


        virtual scalar unscaledEnergy(const scalar r) const = 0;

        virtual const dictionary& dict() const = 0;

        scalar scaledEnergy
        (
            const scalar r
        ) const;

        scalar energyDerivative
        (
            const scalar r,
            const bool scaledEnergyDerivative = true
        ) const;

        const dictionary& pairPotentialProperties() const;

        void output(const fileName& pathName);
        
        void writeEnergyAndForceTables(const fileName& pathName);

        //- Read pairPotential dictionary
//         virtual bool read
//         (
//             const dictionary& pairPotentialProperties
//         ) = 0;

        virtual void write(const fileName& pathName) = 0;
        

        inline const word& name() const;
        inline const List<word>& idList() const;
        inline scalar rMin() const;
        inline scalar dr() const;
        inline scalar rCut() const;
        inline scalar rCutSqr() const;        
        inline bool writeTables() const;        
        inline bool exclusions() const;    
        
        inline autoPtr<exclusionModel>& excludeModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "pairPotentialModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
