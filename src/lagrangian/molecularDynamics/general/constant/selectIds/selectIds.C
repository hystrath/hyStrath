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
    selectIds

Description

\*----------------------------------------------------------------------------*/

#include "selectIds.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
selectIds::selectIds()
:
    molIds_()
{}


selectIds::selectIds
(
    const constantMoleculeProperties& cP,
    const dictionary& dict
)
:
    molIds_()
{
    const List<word> molecules (dict.lookup("molIds"));
    
    if(molecules.size() > 0)
    {
        if
        ( 
            (molecules.size() == 1) &&
            (molecules[0] == "ALL")
        )
        {
            molIds_.setSize(cP.nMolTypes(), -1);
            
            forAll(cP.molIds(), i)
            {
                molIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> moleculesReduced(0);
        
            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);
        
                if(findIndex(moleculesReduced, moleculeName) == -1)
                {
                    moleculesReduced.append(moleculeName);
                }
            }
        
            molIds_.setSize(moleculesReduced.size(), -1);
            molIdNames_.setSize(moleculesReduced.size());
            
            forAll(moleculesReduced, i)
            {
                const word& moleculeName(moleculesReduced[i]);
        
                label molId(findIndex(cP.molIds(), moleculeName));
        
                if(molId == -1)
                {
                    FatalErrorIn
                    (
                        "selectIds::selectIds()"
                    )
                        << "Cannot find id: " << moleculeName << nl << "in dictionary."
                        << exit(FatalError);
                }
        
                molIds_[i] = molId;
                molIdNames_[i] = moleculeName;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "selectIds::selectIds()"
        )
            << "molIds need to be greater than 0 in dictionary."
            << exit(FatalError);
    }
}

selectIds::selectIds
(
    const constantMoleculeProperties& cP,
    const dictionary& dict,
    const word& molIdsHeader
)
:
    molIds_()
{
    const List<word> molecules (dict.lookup(molIdsHeader));

    if(molecules.size() > 0)
    {
        if
        ( 
            (molecules.size() == 1) &&
            (molecules[0] == "ALL")
        )
        {
            molIds_.setSize(cP.nMolTypes(), -1);
            
            forAll(cP.molIds(), i)
            {
                molIds_[i] = i;
            }
        }
        else
        {
            DynamicList<word> moleculesReduced(0);
        
            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);
        
                if(findIndex(moleculesReduced, moleculeName) == -1)
                {
                    moleculesReduced.append(moleculeName);
                }
            }
        
            molIds_.setSize(moleculesReduced.size(), -1);
            molIdNames_.setSize(moleculesReduced.size());
            
            forAll(moleculesReduced, i)
            {
                const word& moleculeName(moleculesReduced[i]);
        
                label molId(findIndex(cP.molIds(), moleculeName));
        
                if(molId == -1)
                {
                    FatalErrorIn
                    (
                        "selectIds::selectIds()"
                    )
                        << "Cannot find id: " << moleculeName << nl << "in dictionary."
                        << exit(FatalError);
                }
        
                molIds_[i] = molId;
                molIdNames_[i] = moleculeName;
                
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "selectIds::selectIds()"
        )
            << "molIds need to be greater than 0 in dictionary."
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

selectIds::~selectIds()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// void selectIds::operator=(const selectIds& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("selectIds::operator=(const selectIds&)")
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

// Ostream& operator<<(Ostream& os, const selectIds& d)
// {
//     os  << d.binWidth_
//         << static_cast<const Map<label>&>(d);
// 
//     // Check state of Ostream
//     os.check
//     (
//         "Ostream& operator<<(Ostream&, "
//         "const selectIds&)"
//     );
// 
//     return os;
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
