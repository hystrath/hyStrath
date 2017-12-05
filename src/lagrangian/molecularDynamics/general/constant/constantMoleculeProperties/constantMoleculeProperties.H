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
    constantMoleculeProperties

Description


SourceFiles
    constantMoleculeProperties.C

\*---------------------------------------------------------------------------*/

#ifndef constantMoleculeProperties_H
#define constantMoleculeProperties_H

#include "Time.H"
#include "polyMesh.H"
#include "reducedUnits.H"
#include "diagTensor.H"
#include "transform.H"
// #include "rotationTensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constantMoleculeProperties Declaration
\*---------------------------------------------------------------------------*/

class constantMoleculeProperties
{
    // Private data

        const polyMesh& mesh_;
        
        const reducedUnits& rU_;
        
        // number of molecules types 
        label N_;
        
        // ordered list of molecule names
        List<word> names_;
        
        // list of total masses of each molecule type
        List<scalar> masses_;
        
        // list of moment of inertia of each molecule type
        List<diagTensor> momentOfInertia_;
        
        // outer list = molecule types, linked to cloud (e.g. polyMoleculeCloud)
        List<word> cloudTypes_;

        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = name
        List<List<word> > siteNames_;

        
        // outer list = molecule types
        // inner list = all sites in molecule
        // entry true/false if site is a pair potential
        List<List<bool> > pairPots_;        
        
        // outer list = molecule types
        // inner list = only pair potential sites
        // entry = link to siteNames_
        List<List<label> > pairSitesLinks_;
        
        // outer list = molecule types
        // inner list = just pair-pot sites in molecule // duplicates allowed
        // entry = name
        List<List<word> > pairPotNames_;

     
        
        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = position of site
        List<List<vector> > siteRefPositions_;

        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = site mass
        List<List<scalar> > siteMasses_;
        
        
        // outer list = molecule types
        // inner list = all sites in molecule 
        // entry = site charge
        List<List<scalar> > siteCharges_;

        // outer list = molecule types
        // inner list = only charged sites
        // entry = link to siteNames_
        List<List<label> > chargeSitesLinks_;        

        // outer list = molecule types
        // inner list = just charge sites in molecule // duplicates allowed
        // entry = name
        List<List<word> > chargeNames_;   
        
        
        label nSites_;
        label nPairPotSites_;
        label nChargeSites_;        
        
        // ordered unqiue lists (A,B,C, ..., N)
       
        // full list of sites (unique list)
        // Outer list = total number of sites in system
        // Entry = name
        List<word> siteIdList_;
        
        // just those that are pair potential (unique list)
        // Outer list = total number of pair pot sites in system
        // Entry = name          // duplicates not allowed      
        List<word> pairPotSiteIdList_;
        
        // just those sites that are charged (unique list)
        // Outer list = total number of non-zero charged sites in system
        // Entry = name        // duplicates not allowed
        List<word> chargeSiteIdList_;


        // linked lists
        
        // link between siteNames_ and siteIdList_
        // outer list = molecule types
        // inner list = all sites in molecule
        // entry = index to siteIdList_
        List<List<label> > siteNames_to_siteIdList_;
        
        
        
public:

    // Constructors

        //- constructor
        constantMoleculeProperties
        (
            const polyMesh& mesh,
            const reducedUnits& rU
        );

    // Destructor

        ~constantMoleculeProperties();


    // Member functions

        bool linearMoleculeTest(const label& idI) const;
        
        const List<diagTensor>& momentOfInertia() const;
        const diagTensor& momentOfInertia(const label& id) const;
        bool linearMolecule(const label& id) const;
        bool pointMolecule(const label& id) const;
        label degreesOfFreedom(const label& id) const;
        
        // no of id types
        const label& nMolTypes() const;
        
        // namelist of all molecule ids 
        const List<word>& molIds() const;
        
        // namelist of cloud types (redundant, unless various clouds exist)
        const List<word>& cloudTypes() const;
        
        
        const List<scalar>& mass() const;
        const scalar& mass(const label& id) const;        
        
        label nSites(const word& idName) const;
        label nSites(const label& id) const;
        
        const List<List<word> >& siteNames() const;        
        const List<word>& siteNames(const word& idName) const;
        const List<word>& siteNames(const label& id) const;

        const List<word>& pairPotNames(const word& idName) const;
        const List<word>& pairPotNames(const label& id) const;
        
        const List<List<label> >& pairSitesLinks() const;
        
        const List<List<vector> >& siteRefPositions() const;
        const List<vector>& siteRefPositions(const word& idName) const;
        const List<vector>& siteRefPositions(const label& id) const;

        const List<scalar>& siteMasses(const word& idName) const;
        const List<scalar>& siteMasses(const label& id) const;
        
        const List<scalar>& siteCharges(const word& idName) const;
        const List<scalar>& siteCharges(const label& id) const;  
        
        const List<List<label> >& chargeSitesLinks() const;
        
        const label& nSiteTypes() const;
        const label& nPairPotTypes() const;
        const label& nChargeTypes() const;
        
        // unique lists of sites
        const List<word>& siteIds() const;
        const List<word>& pairPotSiteIdList() const;
        const List<word>& chargeSiteIdList() const;
        
        
        
        const List<List<label> >& siteNames_to_siteIdList() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
