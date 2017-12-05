/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    potentials

Description
    The first point of contact of storage and access of potentials 
    

SourceFiles
    potentialsI.H
    potentials.C

\*---------------------------------------------------------------------------*/

#ifndef potentials_H
#define potentials_H

#include "polyMesh.H"
#include "IOdictionary.H"
#include "Time.H"
#include "reducedUnits.H"
#include "constantMoleculeProperties.H"
#include "pairPotentials.H"
#include "polyMolecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class potentials Declaration
\*---------------------------------------------------------------------------*/

class potentials
{
    // Private data

        const polyMesh& mesh_;
        
        polyMoleculeCloud& molCloud_;

        const reducedUnits& redUnits_;
        
        const constantMoleculeProperties& cP_;
        
        pairPotentials pairPotentials_;
        
        scalar potentialEnergyLimit_;

        labelList removalOrder_;
        
        // largest rCut of all pair and electrostatic potentials
        scalar rCutMax_; 
        
        
        bool checkPotentialOverlaps_;
        
        // add exclusions functionality here 
        

        // user defined potentials go here and chosen from potentialsDict
        
        
        
        
    // Private Member Functions

        void readPotentialDict();

        
        //- Disallow default bitwise copy construct
        potentials(const potentials&);

        //- Disallow default bitwise assignment
        void operator=(const potentials&);


public:

        friend class polyMolecule;
        
    // Constructors

        //- Construct from standard quantities
        potentials
        (
            const polyMesh& mesh,
            polyMoleculeCloud& molCloud,
            const reducedUnits& rU,
            const constantMoleculeProperties& cP
        );



    //- Destructor
        ~potentials();


    // Member Functions

        // Access

//             inline label nIds() const;

//             inline const List<word>& idList() const;

//             inline const List<word>& siteIdList() const;

//             bool interactSites(polyMolecule* molI, polyMolecule* molJ, siteI, siteJ);
            
//             const pairPotentials& pairPots() const;
            pairPotentials& pairPots();            
//             pairPotentials& pairPotentials();
            
            inline const scalar& potentialEnergyLimit() const;

            inline const labelList& removalOrder() const;

//             inline const pairPotentialList& pairPotentials() const;
            
            const scalar& rCutMax() const;
            
            bool checkPotentialOverlaps();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "potentialsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
