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

Diffuse wall patch interaction for which atoms and atoms ions are fully
catalyzed at the surface to their respective molecules.

\*---------------------------------------------------------------------------*/

#include "dsmcFullyCatalyticDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcFullyCatalyticDiffuseWallPatch, 0);

addToRunTimeSelectionTable
( 
    dsmcPatchBoundary, 
    dsmcFullyCatalyticDiffuseWallPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcFullyCatalyticDiffuseWallPatch::setProperties()
{
    dsmcDiffuseWallPatch::setProperties();
    dsmcCatalyticWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFullyCatalyticDiffuseWallPatch::dsmcFullyCatalyticDiffuseWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseWallPatch(t, mesh, cloud, dict),
    dsmcCatalyticWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFullyCatalyticDiffuseWallPatch::~dsmcFullyCatalyticDiffuseWallPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcFullyCatalyticDiffuseWallPatch::initialConfiguration()
{}


void dsmcFullyCatalyticDiffuseWallPatch::calculateProperties()
{}


void dsmcFullyCatalyticDiffuseWallPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    const label& iD = findIndex(catalysisTypeIds_, p.typeId());
    
    if(iD == -1)
    {       
        //- diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection(p);
    }
    else
    {          
        //- Edit the particle's identity
        p.typeId() = catalysedTypeIds_[iD];
        
        //- diffuse reflection
        dsmcDiffuseWallPatch::performDiffuseReflection(p);

        measurePropertiesAfterControl(p, heatOfReaction_[iD]);
    }
}


void dsmcFullyCatalyticDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcFullyCatalyticDiffuseWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
