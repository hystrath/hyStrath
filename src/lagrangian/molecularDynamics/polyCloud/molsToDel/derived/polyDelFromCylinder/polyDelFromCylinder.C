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

#include "polyDelFromCylinder.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDelFromCylinder, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyDelFromCylinder, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDelFromCylinder::polyDelFromCylinder
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    radius_(readScalar(propsDict_.lookup("radius")))

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

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    bool invert = false;

    if (propsDict_.found("invert"))
    {
        invert = Switch(propsDict_.lookup("invert"));
    }    
    
    if(invert)
    {
        findMolsToDelOutside();
    }
    else
    {
        findMolsToDelInside();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDelFromCylinder::~polyDelFromCylinder()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyDelFromCylinder::findMolsToDelInside()
{
    DynamicList<polyMolecule*> molsToDel;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    label initialSize = molCloud_.size();

    scalar rSEMag = mag(endPoint_ - startPoint_);

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        const vector& rI = mol().position();
        vector rSI = rI - startPoint_;
        scalar centreLineDistance = rSI & unitVector_;

        //- step 1: test polyMolecule is between starting point and end point
        if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
        {
            vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

            //step 2: test polyMolecule is within radial distance of centre-line
            if(mag(pointOnCentreLine-rI) <= radius_)
            {
                label molId = mol().id();

                if(findIndex(molIds_, molId) != -1)
                {
                    polyMolecule* molI = &mol();
                    molsToDel.append(molI);
                }
            }
        }
    }

    //molsToDel.shrink();

    forAll(molsToDel, m)
    {        
        Info <<  molsToDel[m]->position() << endl;
        
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
}

void polyDelFromCylinder::findMolsToDelOutside()
{
    DynamicList<polyMolecule*> molsToDel;

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    label initialSize = molCloud_.size();

    scalar rSEMag = mag(endPoint_ - startPoint_);

    for
    (
        mol = molCloud_.begin();
        mol != molCloud_.end();
        ++mol
    )
    {
        const vector& rI = mol().position();
        vector rSI = rI - startPoint_;
        scalar centreLineDistance = rSI & unitVector_;
        
        bool del = true;
         
        if(findIndex(molIds_, mol().id()) == -1)
        {
            del = false;
        }

 
        //- step 1: test polyMolecule is between starting point and end point
        if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
        {
            vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

            //step 2: test polyMolecule is within radial distance of centre-line
            if(mag(pointOnCentreLine-rI) <= radius_)
            {
                //label molId = mol().id();

                //if(findIndex(molIds_, molId) != -1)
                {
                    del = false;
                }
            }
        }
        
        if(del)
        {
            polyMolecule* molI = &mol();
            molsToDel.append(molI);
        }
    }

    molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
}


} // End namespace Foam

// ************************************************************************* //
