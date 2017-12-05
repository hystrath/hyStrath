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

Abstract class for wall patch catalysis.

\*---------------------------------------------------------------------------*/

#include "dsmcCatalyticWallPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcCatalyticWallPatch, 0);

// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcCatalyticWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));
    
    // set the molecules/atoms to be catalysed typeIds ------------    
    const List<word> inputMolecules 
            (propsDict_.lookup("moleculesToBeCatalysed"));
            
    if(inputMolecules.size() == 0)
    {
        
        FatalErrorIn("dsmcCatalyticWallPatch::setProperties()")
            << "Cannot have zero typeIds being catalysed." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }
    
    DynamicList<word> inputMoleculesReduced(0);

    forAll(inputMolecules, i)
    {
        const word& moleculeName(inputMolecules[i]);

        if(findIndex(inputMoleculesReduced, moleculeName) == -1)
        {
            inputMoleculesReduced.append(moleculeName);
        }
    }

    inputMoleculesReduced.shrink();
    
    //  set the type ids
    
    catalysisTypeIds_.setSize(inputMoleculesReduced.size(), -1); 
    
    forAll(inputMoleculesReduced, i)
    {
        const word& moleculeName(inputMoleculesReduced[i]);
        
        label typeId = findIndex(cloud_.typeIdList(), moleculeName);
        
        // check that input molecules belong to the typeIdList
        if(typeId == -1)
        {
            FatalErrorIn("dsmcCatalyticWallPatch::setProperties()")
                << "Cannot find type id: " << moleculeName << nl 
                << exit(FatalError);
        }
        
        catalysisTypeIds_[i] = typeId;
    }
    
    // set the product molecules/atoms typeIds   
    
    const List<word> outputMolecules (propsDict_.lookup("catalysedMolecules"));
    
    if(outputMolecules.size() != inputMolecules.size())
    {
        
        FatalErrorIn("dsmcCatalyticWallPatch::setProperties()")
            << "catalysedMolecules must be the same size as "
            << "moleculesToBeCatalysed." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }
    
    DynamicList<word> outputMoleculesReduced(0);

    forAll(outputMolecules, i)
    {
        const word& moleculeName(outputMolecules[i]);

        if(findIndex(outputMoleculesReduced, moleculeName) == -1)
        {
            outputMoleculesReduced.append(moleculeName);
        }
    }

    outputMoleculesReduced.shrink();
       
    catalysedTypeIds_.setSize(outputMoleculesReduced.size(), -1);
    
    forAll(outputMoleculesReduced, i)
    {
        const word& moleculeName(outputMoleculesReduced[i]);
        
        label typeId = findIndex(cloud_.typeIdList(), moleculeName);
        
        // check that output molecules belong to the typeIdList
        if(typeId == -1)
        {
            FatalErrorIn("dsmcCatalyticWallPatch::setProperties()")
                << "Cannot find type id: " << moleculeName << nl 
                << exit(FatalError);
        }
        
        catalysedTypeIds_[i] = typeId;
    }
    
    heatOfReaction_.clear();
    
    heatOfReaction_.setSize(catalysisTypeIds_.size(), 0.0);
    
    if(heatOfReaction_.size() == inputMolecules.size())
    {
        // set the heat of reactions
       
        forAll(heatOfReaction_, i)
        {
            heatOfReaction_[i] = readScalar
            (
                propsDict_.subDict("heatOfReaction").lookup(inputMoleculesReduced[i])
            );
        }
    }
    else
    {
        FatalErrorIn("dsmcCatalyticWallPatch::setProperties()")
                << "heatOfReaction list must be same size as "
                << "moleculesToBeCatalysed." << nl 
                << exit(FatalError);
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCatalyticWallPatch::dsmcCatalyticWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    catalysisTypeIds_(),
    catalysedTypeIds_(),
    heatOfReaction_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCatalyticWallPatch::~dsmcCatalyticWallPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcCatalyticWallPatch::initialConfiguration()
{}


void dsmcCatalyticWallPatch::calculateProperties()
{}


void dsmcCatalyticWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcCatalyticWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
