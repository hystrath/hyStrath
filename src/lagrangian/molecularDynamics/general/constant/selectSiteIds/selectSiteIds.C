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
    selectSiteIds

Description

\*----------------------------------------------------------------------------*/

#include "selectSiteIds.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
selectSiteIds::selectSiteIds()
:
    siteIds_()
{}



selectSiteIds::selectSiteIds
(
    const constantMoleculeProperties& cP,
    const dictionary& dict
)
:
    siteIds_()
{
    const List<word> sites (dict.lookup("siteIds"));

    if(sites.size() > 0)
    {
        if
        ( 
            (sites.size() == 1) &&
            (sites[0] == "ALL")
        )
        {
            siteIds_.setSize(cP.nSiteTypes(), -1);
            
            forAll(cP.siteIds(), i)
            {
                siteIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> sitesReduced(0);
        
            forAll(sites, i)
            {
                const word& siteName(sites[i]);
        
                if(findIndex(sitesReduced, siteName) == -1)
                {
                    sitesReduced.append(siteName);
                }
            }
        
            siteIds_.setSize(sitesReduced.size(), -1);
            siteIdNames_.setSize(sitesReduced.size());
            
            forAll(sitesReduced, i)
            {
                const word& siteName(sitesReduced[i]);
        
                label siteId(findIndex(cP.siteIds(), siteName));
        
                if(siteId == -1)
                {
                    FatalErrorIn
                    (
                        "selectSiteIds::selectSiteIds()"
                    )
                        << "Cannot find siteId: " << siteName
                        << " in siteIds dictionary = " << sites
                        << nl << exit(FatalError);
                }
        
                siteIds_[i] = siteId;
                siteIdNames_[i] = siteName;
            }
        }
    }
//     else
//     {
//         FatalErrorIn
//         (
//             "selectSiteIds::selectSiteIds()"
//         )
//             << "siteIds need to be greater than 0 in dictionary."
//             << exit(FatalError);
//     }
}

selectSiteIds::selectSiteIds
(
    const constantMoleculeProperties& cP,
    const dictionary& dict,
    const word& siteIdsHeader
)
:
    siteIds_()
{
    const List<word> sites (dict.lookup(siteIdsHeader));

    if(sites.size() > 0)
    {
        if
        ( 
            (sites.size() == 1) &&
            (sites[0] == "ALL")
        )
        {
            siteIds_.setSize(cP.nSiteTypes(), -1);
            
            forAll(cP.siteIds(), i)
            {
                siteIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> sitesReduced(0);
        
            forAll(sites, i)
            {
                const word& siteName(sites[i]);
        
                if(findIndex(sitesReduced, siteName) == -1)
                {
                    sitesReduced.append(siteName);
                }
            }
        
            siteIds_.setSize(sitesReduced.size(), -1);
            siteIdNames_.setSize(sitesReduced.size());
            
            forAll(sitesReduced, i)
            {
                const word& siteName(sitesReduced[i]);
        
                label siteId(findIndex(cP.siteIds(), siteName));
        
                if(siteId == -1)
                {
                    FatalErrorIn
                    (
                        "selectSiteIds::selectSiteIds()"
                    )
                        << "Cannot find siteId: " << siteName
                        << " in siteIds dictionary = " << sites
                        << nl << exit(FatalError);
                }
        
                siteIds_[i] = siteId;
                siteIdNames_[i] = siteName;                
            }
        }
    }
//     else
//     {
//         FatalErrorIn
//         (
//             "selectSiteIds::selectSiteIds()"
//         )
//             << "siteIds need to be greater than 0 in dictionary."
//             << exit(FatalError);
//     }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

selectSiteIds::~selectSiteIds()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void selectSiteIds::operator=(const selectSiteIds& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("selectSiteIds::operator=(const selectSiteIds&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// 
//     Map<label>::operator=(rhs);
// 
//     binWidth_ = rhs.binWidth();
// }


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Ostream& operator<<(Ostream& os, const selectSiteIds& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const selectSiteIds&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
