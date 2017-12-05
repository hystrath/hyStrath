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
// void selectIds::setRadius()
// {
//     for(label i = 0; i < noOfBins_; i++)
//     {
//        radius_[i] = (0.5 + scalar(i)) * binWidth();
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
selectIds::selectIds()
:
    molIds_()
{}



selectIds::selectIds
(
    const potential& pot,
    const dictionary& dict
)
:
    molIds_()
{
    const List<word> molecules (dict.lookup("molIds"));

    if(molecules.size() > 0)
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
    
        forAll(moleculesReduced, i)
        {
            const word& moleculeName(moleculesReduced[i]);
    
            label molId(findIndex(pot.idList(), moleculeName));
    
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
        }
    }
    else
    {
        molIds_.setSize(pot.idList().size(), -1);

        forAll(molIds_, i)
        {
            molIds_[i] = i;
        }
    }
}

selectIds::selectIds
(
    const potential& pot,
    const dictionary& dict,
    const word& molIdsName
)
:
    molIds_()
{
    const List<word> molecules (dict.lookup(molIdsName));

    if(molecules.size() > 0)
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
    
        forAll(moleculesReduced, i)
        {
            const word& moleculeName(moleculesReduced[i]);
    
            label molId(findIndex(pot.idList(), moleculeName));
    
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
        }
    }
    else
    {
        molIds_.setSize(pot.idList().size(), -1);

        forAll(molIds_, i)
        {
            molIds_[i] = i;
        }
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
