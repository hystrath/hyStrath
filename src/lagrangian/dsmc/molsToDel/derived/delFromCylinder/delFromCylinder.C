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

#include "delFromCylinder.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(delFromCylinder, 0);

addToRunTimeSelectionTable(molsToDeleteModel, delFromCylinder, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
delFromCylinder::delFromCylinder
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    molsToDeleteModel(cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    radius_(readScalar(propsDict_.lookup("radius"))),
    typeIds_()
{

    // check if start point is in the mesh
   
    if(mesh_.findCell(startPoint_) == -1)
    {
        Info<< "WARNING: starting point " << startPoint_ 
            << " is selected outside the mesh."
            << endl;
    }

    if(mesh_.findCell(endPoint_) == -1)
    {
        Info<< "WARNING: end point " << endPoint_ 
            << " is selected outside the mesh."
            << endl;
    }


    // standard to reading typeIds ------------ 
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

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcInflowPatch::dsmcInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"molsToDeleteDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------


    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

delFromCylinder::~delFromCylinder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void delFromCylinder::findMolsToDel()
{
    DynamicList<dsmcParcel*> molsToDel;

    label initialSize = cloud_.size();

    scalar rSEMag = mag(endPoint_ - startPoint_);

    forAllIter(dsmcCloud, cloud_, p)
    {
//         dsmcParcel& p = iter();

        const vector& rI = p().position();
        vector rSI = rI - startPoint_;
        scalar centreLineDistance = rSI & unitVector_;

        //- step 1: test dsmcParcel is between starting point and end point
        if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
        {
            vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

            //step 2: test dsmcParcel is within radial distance of centre-line
            if(mag(pointOnCentreLine-rI) <= radius_)
            {
                if(findIndex(typeIds_, p().typeId()) != -1)
                {
                    dsmcParcel* pI = &p();
                    molsToDel.append(pI);
                }
            }
        }
    }

    molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial molecules: " <<  initialSize 
        << ", molecules kept: " <<  molsKept
        << ", molecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
//     cloud_.rebuildCellOccupancy();
}


} // End namespace Foam

// ************************************************************************* //
