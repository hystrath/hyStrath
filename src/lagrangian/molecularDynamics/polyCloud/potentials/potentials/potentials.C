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

\*---------------------------------------------------------------------------*/

#include "potentials.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void potentials::readPotentialDict()
{
    Info<< nl <<  "Reading potentials dictionary:" << endl;

    IOdictionary potentialsDict
    (
        IOobject
        (
            "potentialDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    potentialEnergyLimit_ = readScalar
    (
        potentialsDict.lookup("potentialEnergyLimit")
    );


    if(redUnits_.runReducedUnits())
    {
        potentialEnergyLimit_ /= redUnits_.refEnergy();
    }

    // removal order in the case of overlaps
    if (potentialsDict.found("removalOrder"))
    {
        List<word> remOrd = potentialsDict.lookup("removalOrder");

        DynamicList<label> removalOrder(0);

        forAll(remOrd, rO)
        {
            const label id = findIndex(cP_.molIds(), remOrd[rO]);

            if( (id != -1) && (findIndex(removalOrder, id) == -1) )
            {
                removalOrder.append(id);
            }
        }

        // fill in remaining spots arbitrarily
        forAll(cP_.molIds(), i)
        {
            if(findIndex(removalOrder, i) == -1)
            {
                removalOrder.append(i);
            }
        }

        removalOrder_.transfer(removalOrder);
        
        Info << "setup removalOrder = " << removalOrder_ << endl;
    }
    else
    {
        FatalErrorIn("potentials::readPotentialDict()")
            << "removalOrder list not found in system/potentialDict"
            << abort(FatalError);        
        
    }
    
    if (potentialsDict.found("checkOverlaps"))
    {
        checkPotentialOverlaps_ = Switch(potentialsDict.lookup("checkOverlaps"));  
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct to run MD simulation (idList is read in from constant dir)
potentials::potentials
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const reducedUnits& rU,
    const constantMoleculeProperties& cP
)
:
    mesh_(mesh),
    molCloud_(molCloud),
    redUnits_(rU),
    cP_(cP),
    pairPotentials_(mesh, molCloud, cP, rU),
    rCutMax_(pairPotentials_.maxRCut()),
    checkPotentialOverlaps_(true)
{
    readPotentialDict();
    
    // set exclusions
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

potentials::~potentials()
{}



// bool potentials::interactMolecules(polyMolecule* molI, polyMolecule* molJ)
// {
//     if(!molI->frozen() || !molJ->frozen())
//     {
//         return false;
//     }
//    
//     if()
// }
//     if(!molI->frozen() || !molJ->frozen())


// const pairPotentials& potentials::pairPots() const
// {
//     return pairPotentials_;
// }


pairPotentials& potentials::pairPots()
{
    return pairPotentials_;
}

const scalar& potentials::rCutMax() const
{
    return rCutMax_;
}


bool potentials::checkPotentialOverlaps()
{
    return checkPotentialOverlaps_;
}

// pairPotentials& potentials::pairPots()
// {
//     return pairPotentials_;
// }

// Foam::pairPotentials& Foam::potentials::pairPotentials()
// {
//     return pairPotentials_;
// }

} // End namespace Foam

// ************************************************************************* //
