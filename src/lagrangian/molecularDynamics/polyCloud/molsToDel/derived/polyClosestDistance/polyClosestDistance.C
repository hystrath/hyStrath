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

#include "polyClosestDistance.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyClosestDistance, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyClosestDistance, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyClosestDistance::polyClosestDistance
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    distance_(readScalar(propsDict_.lookup("distance")))    
{
 
    {
        const word molIdName = propsDict_.lookup("molIdToDelete");

        label molId(findIndex(molCloud_.cP().molIds(), molIdName));

        if(molId == -1)
        {
            FatalErrorIn
            (
                "polyClosestDistance::polyClosestDistance()"
            )
                << "Cannot find id: " << molIdName << nl << "in dictionary."
                << exit(FatalError);
        }
        
        molToDeleteId_ = molId;
    }

    {
        const word molIdName = propsDict_.lookup("molIdReference");

        label molId(findIndex(molCloud_.cP().molIds(), molIdName));

        if(molId == -1)
        {
            FatalErrorIn
            (
                "polyClosestDistance::polyClosestDistance()"
            )
                << "Cannot find id: " << molIdName << nl << "in dictionary."
                << exit(FatalError);
        }
        
        refMolId_ = molId;
    }    
    
    findMolsToDel();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyClosestDistance::~polyClosestDistance()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyClosestDistance::findMolsToDel()
{
    DynamicList<polyMolecule*> molsToDel;
    DynamicList<label> trackingNumbers;
    
    label initialSize = molCloud_.size();
    
    
    {    
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());
        IDLList<polyMolecule>::iterator molJ(molCloud_.begin());
        
        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
           
            for
            (
                molJ = molCloud_.begin();
                molJ != molCloud_.end();
                ++molJ
            )        
            {
                
                if
                (
                    (molI().id() == refMolId_) &&
                    (molJ().id() == molToDeleteId_) &&
                    (mag(molJ().position() - molI().position()) < distance_)
                )
                {    

                    if(findIndex(trackingNumbers, molJ().trackingNumber()) == -1)
                    {
                        polyMolecule* mol = &molJ();
                        molsToDel.append(mol);
                        trackingNumbers.append(molJ().trackingNumber());
                    }
                }
            }
        }
    }
    
    //molsToDel.shrink();

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
    molCloud_.prepareInteractions();
}

void polyClosestDistance::checkBoundBox
(
    boundBox& b,
    const vector& startPoint,
    const vector& endPoint
)
{
    vector& vMin = b.min();
    vector& vMax = b.max();

    if(startPoint.x() < endPoint.x())
    {
        vMin.x() = startPoint.x();
        vMax.x() = endPoint.x();
    }
    else
    {
        vMin.x() = endPoint.x();
        vMax.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        vMin.y() = startPoint.y();
        vMax.y() = endPoint.y();
    }
    else
    {
        vMin.y() = endPoint.y();
        vMax.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        vMin.z() = startPoint.z();
        vMax.z() = endPoint.z();
    }
    else
    {
        vMin.z() = endPoint.z();
        vMax.z() = startPoint.z();
    }
}
} // End namespace Foam

// ************************************************************************* //
