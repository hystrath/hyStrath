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

Description

\*---------------------------------------------------------------------------*/

#include "pairPotentials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// constructor
pairPotentials::pairPotentials
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const constantMoleculeProperties& cP,
    const reducedUnits& redUnits
)
:
    mesh_(mesh),
    cP_(cP),
    potentialsDict_
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    pairPotentialsList_(potentialsDict_.lookup("pairs")),
    pairNames_(pairPotentialsList_.size()),
    pairIds_(pairPotentialsList_.size()),
    pairPotentials_(pairPotentialsList_.size()),
    exclusions_(pairPotentialsList_.size())

{
   
    if(pairPotentialsList_.size() > 0 )
    {
        Info << nl << "Creating potentials: " << nl << endl;
    
        forAll(pairPotentialsList_, i)
        {
            const entry& potI = pairPotentialsList_[i];
            const dictionary& potIDict = potI.dict();
            const word& headerName=potI.keyword();
            Info << "headerName = " << headerName << endl;
            pairPotentials_[i] = autoPtr<pairPotentialModel>
            (
                pairPotentialModel::New(mesh, molCloud, redUnits, headerName, potIDict)
            );
    
            pairNames_[i] = pairPotentials_[i]->type();
            pairIds_[i] = i;
            
            exclusions_[i] = pairPotentials_[i]->exclusions();
        }
        
        testPairPotentials();        
    }
    else
    {
         FatalErrorIn("pairPotentials.C") << nl
                << " Something went wrong. You have no potentials in your system/potentialDict"
                << nl << nl << "Check that the format is as follows: " 
                << " pairs " << nl
                << " ( " << nl
                << "    enter your pair potentials here "
                << " ); " 
                << nl << abort(FatalError);       
    }
    
    // making directory
    pathName_ = mesh_.time().path()/"potentials";

    if(isDir(pathName_) )
    {
        rmDir(pathName_);
    }     
    
    mkDir(pathName_);    
   
    // electrostatic potential 
    
    word headerName = "electrostatic";
    
    dictionary& dict = potentialsDict_.subDict(headerName);
    
    electrostaticPotential_ = autoPtr<pairPotentialModel>
    (
        pairPotentialModel::New(mesh, molCloud, redUnits, headerName, dict)
    );

    Info << "exclusions = " << exclusions_ << endl;
    
    rCut_ = maxRCut();
//     rCutSqr_ = rCut_*rCut_;
    
    
    // write out potentials
    forAll(pairPotentials_, i)
    {
        pairPotentials_[i]->output(pathName_);
    }
    
    electrostaticPotential_->output(pathName_);

}


pairPotentials::~pairPotentials()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void pairPotentials::testPairPotentials()
{
    Info << nl << "Testing pair potentials..." << nl << endl;

    // check that all pair sites have pair potential combinations defined in potentialDict
    // and create the linked lists 

    const List<word>& pairSites = cP_.pairPotSiteIdList();
    
    Info << "cP_.pairPotSiteIdList() " << cP_.pairPotSiteIdList() << endl;
    
    // imp: first to set size of lists
    pairPotIdList_to_pairPotentials_.setSize(pairSites.size());
    
    forAll(pairSites, i)
    {
        pairPotIdList_to_pairPotentials_[i].setSize(pairSites.size(), -1); 
    }
    
    forAll(pairSites, i)
    {
        forAll(pairSites, j)
        {
            bool pairsFound = false;
            
            forAll(pairPotentials_, k)
            {
                const List<word>& ids = pairPotentials_[k]->idList();
                
                if
                (
                    ( (ids[0] == pairSites[i]) && (ids[1] == pairSites[j]) ) ||
                    ( (ids[1] == pairSites[i]) && (ids[0] == pairSites[j]) ) 
                )
                {
                    pairPotIdList_to_pairPotentials_[i][j] = k;
                    pairPotIdList_to_pairPotentials_[j][i] = k;
                    pairsFound = true;
                }
            }
            
            if(!pairsFound)
            {
                FatalErrorIn("pairPotentials.C") << nl
                        << " You did not specify the following pair combination => "
                        << pairSites[i] << "-" << pairSites[j]
                        << nl << " in the system/potentialDict."
                        << nl << abort(FatalError);                 
                
            }
        }
    }
    
    // create linked lists for sites
    
    const List<word>& sites = cP_.siteIds();
    
    siteIdList_to_pairPotentials_.setSize(sites.size());
    
    forAll(sites, i)
    {
        siteIdList_to_pairPotentials_[i].setSize(sites.size(), -1); 
    }
    
    forAll(sites, i)
    {
        forAll(sites, j)
        {
            forAll(pairPotentials_, k)
            {
                const List<word>& ids = pairPotentials_[k]->idList();
                
                if
                (
                    ( (ids[0] == sites[i]) && (ids[1] == sites[j]) ) ||
                    ( (ids[1] == sites[i]) && (ids[0] == sites[j]) ) 
                )
                {
                    siteIdList_to_pairPotentials_[i][j] = k;
                    siteIdList_to_pairPotentials_[j][i] = k;
                }
            }
        }
    }
    
    Info << "siteIdList_to_pairPotentials_ = " << siteIdList_to_pairPotentials_ << endl;
}




scalar pairPotentials::maxRCut()
{
    scalar rCut = 0.0;
    
    forAll(pairPotentials_, i)
    {
        if(pairPotentials_[i]->rCut() > rCut)
        {
            rCut = pairPotentials_[i]->rCut();
        }
    }
    
    if(electrostaticPotential_->rCut() > rCut)
    {
        rCut = electrostaticPotential_->rCut();
    }    

    Info << nl  << "rCut = " << rCut << nl << endl;
    
    return rCut;
}        

scalar pairPotentials::force
(
    const label k,
    const scalar r
) const
{
    return pairPotentials_[k]->force(r);
}

scalar pairPotentials::energy
(
    const label k,
    const scalar r
) const
{
    return pairPotentials_[k]->energy(r);
}


bool pairPotentials::rCutSqr
(
    const label k,
    const scalar rIJSqr
) const
{
    if(rIJSqr < pairPotentials_[k]->rCutSqr())
    {
        return true;
    }
    else
    {
        return false;
    }
}


scalar pairPotentials::rMin
(
    const label& k
) 
{
     return pairPotentials_[k]->rMin();
}

bool pairPotentials::excludeSites
(
    polyMolecule* molI,
    polyMolecule* molJ,
    const label& sI,
    const label& sJ,
    const label& k
) 
{
    if(!exclusions_[k])
    {
        return false;
    }
    else
    {
        return pairPotentials_[k]->excludeModel()->excludeSites(molI, molJ, sI, sJ);
    }
}

void pairPotentials::initialiseExclusionModels()
{
    forAll(pairPotentials_, k)
    {
        if(exclusions_[k])
        {
            pairPotentials_[k]->excludeModel()->initialiseExclusions();
        }
    }
}

// const List< List<label> >& pairPotentials::pairPotIdList_to_pairPotentials() const
// {
//     return pairPotIdList_to_pairPotentials_;
// }



// label pairPotentials::pairPotentialIndex
// (
//     const label idI,
//     const label idJ,
//     const label sI,
//     const label sJ
// ) const
// {
//     
//     
// }

// supply molId I, molId J and i and j pot site indexes (see format of pairPotNames_ in cP_)
//make inline
// label pairPotentials::getIndexFromPairPotentialSites
// (
//     const label idI,
//     const label idJ,
//     const label i,
//     const label j
// ) const
// {
//     label idsI=cP_.pairPotNames_to_pairPotSitesList()[idI][i];
//     label idsJ=cP_.pairPotNames_to_pairPotSitesList()[idJ][j];
//     return pairPotIdList_to_pairPotentials_[idsI][idsJ];
// }

// supply molId I, molId J and sI and sJ site indexes (see format of siteNames_ in cP_)
//make inline

label pairPotentials::pairPotentialIndex
(
    const label idI,
    const label idJ,
    const label sI,
    const label sJ
) const
{
    label idsI=cP_.siteNames_to_siteIdList()[idI][sI];
    label idsJ=cP_.siteNames_to_siteIdList()[idJ][sJ];
    return siteIdList_to_pairPotentials_[idsI][idsJ];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
