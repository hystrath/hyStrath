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
    pairPotentials

Description
    Storage and handling of all pair potentials (including electrostatics)

\*---------------------------------------------------------------------------*/

#ifndef pairPotentials_H
#define pairPotentials_H

#include "vector.H"
#include "List.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "polyMesh.H"
#include "pairPotentialModel.H"
#include "reducedUnits.H"
#include "constantMoleculeProperties.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairPotentials Declaration
\*---------------------------------------------------------------------------*/

class pairPotentials
{

    // Private data

        const polyMesh& mesh_;

        const constantMoleculeProperties& cP_;
        
        IOdictionary potentialsDict_;

        fileName pathName_;
        
        // maximum rCut
        scalar rCut_;
        
       
//         scalar rCutSqr_;
        
        // these are just pair potentials
        
        PtrList<entry> pairPotentialsList_;
        List<word> pairNames_;
        List<label> pairIds_;
        List< autoPtr<pairPotentialModel> > pairPotentials_;

        
        // linking pairPotIdList_ with pairPotentials_
        // outer list size of cP_.pairPotSiteIdList()
        // inner list size of cP_.pairPotSiteIdList()
        // store index of pairPotentials_
        List<List<label> > pairPotIdList_to_pairPotentials_;

        // linking siteIdList_ with pairPotentials_
        // outer list size of cP_.siteIds()
        // inner list size of cP_.siteIds()
        // store index of pairPotentials_        
        List<List<label> > siteIdList_to_pairPotentials_;
        
        
        // this is a fixed electrostatic short range potential applied to all
        // charged sites 
        autoPtr<pairPotentialModel> electrostaticPotential_;
        
        
        
        List<bool> exclusions_;
        
       


        
        
        
        // private member functions 
        
            void testPairPotentials();

//             void setExclusions();
        
public:

    // Constructors
    
    
        //- Standard constructor
        pairPotentials
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const constantMoleculeProperties& cP,
            const reducedUnits& redUnits
        );


        ~pairPotentials();

    // Member Functions

        scalar maxRCut();
        
        inline List< autoPtr<pairPotentialModel> >& pairPots();
        inline const List< autoPtr<pairPotentialModel> >& pairPots() const;
        inline const List<word>& pairNames() const;
        inline const List<label>& pairIds() const;
        inline autoPtr<pairPotentialModel>& electrostaticPot();
        inline const autoPtr<pairPotentialModel>& electrostaticPot() const;        
        
//         label getIndexFromPairPotentialSites
//         (
//             const label idI,
//             const label idj,
//             const label i,
//             const label j
//         ) const; 
//         
//         label getIndexFromSites
//         (
//             const label idI,
//             const label idJ,
//             const label sI,
//             const label sJ
//         ) const;  
        
        label pairPotentialIndex
        (
            const label idI,
            const label idJ,
            const label sI,
            const label sJ
        ) const;

        scalar force
        (
            const label k,
            const scalar r
        ) const;
        
        scalar energy
        (
            const label k,
            const scalar r
        ) const;

        bool rCutSqr
        (
            const label k,
            const scalar r
        ) const;         

        scalar rMin
        (
            const label& k
        );
        
        bool excludeSites
        (
            polyMolecule* molI,
            polyMolecule* molJ,
            const label& sI,
            const label& sJ,
            const label& k
        );
        
        
//         const List< List<label> >& pairPotIdList_to_pairPotentials() const;
        
        
        void initialiseExclusionModels();
        
//         
    // Ostream Operator

//         friend Ostream& operator<<(Ostream&, const pairPotentials&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pairPotentialsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
