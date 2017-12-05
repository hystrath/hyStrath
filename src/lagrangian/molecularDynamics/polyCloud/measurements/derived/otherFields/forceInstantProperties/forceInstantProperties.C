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

#include "forceInstantProperties.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(forceInstantProperties, 0);

addToRunTimeSelectionTable(polyField, forceInstantProperties, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void forceInstantProperties::setBoundBoxes()
// {
//  
//     PtrList<entry> boxList(propsDict_.lookup("boxes"));
// 
//     boxes_.setSize(boxList.size());
// 
//     forAll(boxList, b)
//     {
//         const entry& boxI = boxList[b];
//         const dictionary& dict = boxI.dict();
// 
//         vector startPoint = dict.lookup("startPoint");
//         vector endPoint = dict.lookup("endPoint");
//         boxes_[b].resetBoundedBox(startPoint, endPoint);
//     }
// }



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
forceInstantProperties::forceInstantProperties
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName"))
//     boxes_(),
//     molIds_()
{
        // build bound boxes
//     setBoundBoxes();
    
    measureInterForcesSites_ = true;
    
    {
        // choose molecule ids to sample
        molIdsWall_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsWall"
            
        );

        molIdsWall_ = ids.molIds();
    }
    
    {
        // choose molecule ids to sample
        molIdsFluid_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsFluid"
            
        );

        molIdsFluid_ = ids.molIds();
    }
    
    forceField_.clear();
    pairsField_.clear();
    force_ = vector::zero;
    pairs_ = 0.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

forceInstantProperties::~forceInstantProperties()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forceInstantProperties::createField()
{
 
}

void forceInstantProperties::calculateField()
{
    // - parallel processing
    if(Pstream::parRun())
    {
        reduce(force_, sumOp<vector>());
        reduce(pairs_, sumOp<scalar>());
    }
    
    forceField_.append(force_);
    pairsField_.append(pairs_);
    
    force_ = vector::zero;
    pairs_ = 0.0;
}


void forceInstantProperties::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            forceField_.shrink();

            scalarField timeField (forceField_.size(), 0.0);
            vectorField force (forceField_.size(), vector::zero);
            scalarField pairs (pairsField_.size(), 0.0);            
            
            force.transfer(forceField_);
            forceField_.clear();

            pairs.transfer(pairsField_);
            pairsField_.clear();
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "force_instant_"+fieldName_+"_force_SI.xyz",
                timeField*molCloud_.redUnits().refTime(),
                force*molCloud_.redUnits().refForce(),
                true
            );
            
            writeTimeData
            (
                casePath_,
                "force_instant_"+fieldName_+"_pairs.xyz",
                timeField,
                pairs,
                true
            );            
        }
    }
}
void forceInstantProperties::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void forceInstantProperties::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{
    label idI = molI->id();
    label idJ = molJ->id();
    
    label idsI=molCloud_.cP().pairPotNamesToPairPotSitesList()[idI][i];
    label idsJ=molCloud_.cP().pairPotNamesToPairPotSitesList()[idJ][j];  


    label idIW = findIndex(molIdsWall_, molI->id());
    label idJW = findIndex(molIdsWall_, molJ->id());

    label idIF = findIndex(molIdsFluid_, molI->id());
    label idJF = findIndex(molIdsFluid_, molJ->id());
    
//     idIW & idJW = no
//     idIW & idIF = no
//     
//     idIW & idJF = yes
//     idJW & idIF = yes
//     
//     idIF & idJF = no
//     idJF & idJW = no
    
    
    // Assume FLUID-to-WALL forces
    
    if
    (
        ((idIW != -1) && (idJF != -1)) 
    )
    {
        vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
        scalar rsIsJMag = mag(rsIsJ);
        vector force = (rsIsJ/rsIsJMag) * molCloud_.pot().pairPots().force(idsI, idsJ, rsIsJMag);

        
        if(molI->referred() || molJ->referred())
        {
            pairs_ += 0.5;
            force_ += -force*0.5;// This is negative because the wall is particle J
        }
        else
        {
            force_ += -force;// This is negative because the wall is particle J
            pairs_ += 1.0;           
        }        
    }
    else if((idJW != -1) && (idIF != -1))
    {
        vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
        scalar rsIsJMag = mag(rsIsJ);
        vector force = (rsIsJ/rsIsJMag) * molCloud_.pot().pairPots().force(idsI, idsJ, rsIsJMag);
        
        if(molI->referred() || molJ->referred())
        {
            pairs_ += 0.5;
            force_ += force*0.5; 
        }
        else
        {
            force_ += force; 
            pairs_ += 1.0;           
        }
    }    
    
}

const propertyField& forceInstantProperties::fields() const
{
    return fields_;
}

// void forceInstantProperties::updateProperties(const dictionary& newDict)
// {
//     //- the main properties should be updated first
//     updateBasicFieldProperties(newDict);

// }

} // End namespace Foam

// ************************************************************************* //
