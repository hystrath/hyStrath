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

#include "pdAdsorbingWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdAdsorbingWallPatch, 0);

addToRunTimeSelectionTable(pdPatchBoundary, pdAdsorbingWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdAdsorbingWallPatch::pdAdsorbingWallPatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    allSpecies_(false),
    typeElec_(),
    typeIds_()
{
    measurePropertiesAtWall_ = true;
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdAdsorbingWallPatch::~pdAdsorbingWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void pdAdsorbingWallPatch::initialConfiguration()
{}

void pdAdsorbingWallPatch::calculateProperties()
{
    if(typeElec_ == 1) //insulated
    {
        //Info << "Insulated material not implimented!" << endl;
    }
    else if(typeElec_ == 2) // conductor
    {
        scalar totalQ = 0.0;

        pdEmFields& eM = cloud_.emFields();

        forAll(eM.wallQ_.boundaryField()[patchId_],fI)
        {
            totalQ += eM.wallQ_.boundaryField()[patchId_][fI];
        }

        if(Pstream::parRun())
        {
            reduce(totalQ, sumOp<scalar>());  //total charge on all processors
        }

        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        forAll(eM.wallQ_.boundaryField()[patchId_],fI)
        {
            eM.wallQ_.boundaryFieldRef()[patchId_][fI] = totalQ*mag(patch.faceAreas()[fI])/totalPatchSurfaceArea_;
        }
    }
    else
    {
        pdEmFields& eM = cloud_.emFields();
        eM.wallQ_.boundaryFieldRef()[patchId_] = 0;
    }

}

void pdAdsorbingWallPatch::controlParticle(pdParcel& p, pdParcel::trackingData& td)
{

        measurePropertiesBeforeControl(p);

        td.keepParticle = false;

}

void pdAdsorbingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void pdAdsorbingWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}



void pdAdsorbingWallPatch::setProperties()
{
    typeElec_ = readScalar(propsDict_.lookup("typeElec"));

    /*if(typeElec_ == 1)
    {
        Info << "   Warning: " << patchName() << " has been set as insulated (typeElec = 2). Remember to select the appropriate fv boundary scheme." << endl;
    }else if(typeElec_ == 2)
    {
        Info << "   Warning: " << patchName() << " has been set as conducting (typeElec = 2). Remember to select the appropriate fv boundary scheme." << endl;
    }
    else
    {
        Info << "   Warning: " << patchName() << " has been set as grounded. Remember to select the appropriate fv boundary scheme or if this is incorrect, change typeElec (1. insulted, 2. conducting)."  << endl;
    }*/

    if (propsDict_.found("allSpecies"))
    {
        allSpecies_ = Switch(propsDict_.lookup("allSpecies"));
    }

    if(!allSpecies_)
    {
        const List<word> molecules (propsDict_.lookup("typeIds"));

        DynamicList<word> moleculesReduced(0);

        forAll(molecules, i)
        {
            const word& moleculeName(molecules[i]);

            if(findIndex(moleculesReduced, moleculeName) == -1)
            {
                moleculesReduced.append(moleculeName);
            }
        }

        moleculesReduced.shrink();

        typeIds_.clear();

        typeIds_.setSize(moleculesReduced.size(), -1);

        forAll(moleculesReduced, i)
        {
            const word& moleculeName(moleculesReduced[i]);

            label typeId(findIndex(cloud_.typeIdList(), moleculeName));

            if(typeId == -1)
            {
                FatalErrorIn("pdAdsorbingWallPatch::pdAdsorbingWallPatch()")
                    << "Cannot find typeId: " << moleculeName << nl << "in: "
                    << mesh_.time().system()/"boundariesDict"
                    << exit(FatalError);
            }

            typeIds_[i] = typeId;
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
