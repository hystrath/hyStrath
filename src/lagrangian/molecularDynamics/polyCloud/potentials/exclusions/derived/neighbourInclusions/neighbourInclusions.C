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
    neighbourInclusions

Description

\*----------------------------------------------------------------------------*/

#include "neighbourInclusions.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(neighbourInclusions, 0);

addToRunTimeSelectionTable(exclusionModel, neighbourInclusions, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
neighbourInclusions::neighbourInclusions
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    exclusionModel(mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
}


void neighbourInclusions::initialiseExclusions()
{
    // read in tracking numbers 
    List<label> molPointsA = List<label>(propsDict_.lookup("trackingNumbersA"));
    List<label> molPointsB = List<label>(propsDict_.lookup("trackingNumbersB"));
    
    if(molPointsA.size() != molPointsB.size())
    {
        FatalErrorIn("neighbourInclusions::neighbourInclusions()")
            << "size of trackingNumbersA = " << molPointsA
            << " is not the same as trackingNumbersB = : " << molPointsA
            << nl << "in system/potentialDict"
            << exit(FatalError);        
    }
    
    tNsA_.setSize(molPointsA.size());
    tNsB_.setSize(molPointsB.size());
    
    forAll(tNsA_, i)
    {
        tNsA_[i]=molPointsA[i];
        tNsB_[i]=molPointsB[i];
    }
    
    label N=molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    Info << "N = " << N << endl;
    
    fullTNs_.setSize(N);
    
    forAll(tNsA_, i)
    {
        label tNA = tNsA_[i];
        label tNB = tNsB_[i];        
        
        fullTNs_[tNA].append(tNB);
        fullTNs_[tNB].append(tNA);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

neighbourInclusions::~neighbourInclusions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool neighbourInclusions::excludeMolecules
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
    return true;
}

bool neighbourInclusions::excludeSites
(
    polyMolecule* molI,
    polyMolecule* molJ,
    const label& siteI,
    const label& siteJ
)
{
    if(molI->trackingNumber() < fullTNs_.size())
    {
        if(findIndex(fullTNs_[molI->trackingNumber()], molJ->trackingNumber()) == -1)
        {
            return true;
        }
        else
        {
            forAll(fullTNs_[molI->trackingNumber()], i)
            {
                if(fullTNs_[molI->trackingNumber()][i] == molJ->trackingNumber())
                {
                    return false;
                }
            }
            
            return true;
        }
    }
    else
    {
        return true;
    }
}



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
