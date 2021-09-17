/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    newCellZone

Description

\*----------------------------------------------------------------------------*/

#include "newCellZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOdictionary and mesh
newCellZone::newCellZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    cells_()
{
    
    const word zoneName = dict.lookup("zoneName");
    regionName_ = zoneName;
    
    Info << " -> " << zoneName << endl;
        
    const cellZoneMesh& cellZones = mesh.cellZones();
    
    const label& regionId = cellZones.findZoneID(zoneName);
    
    if (regionId != -1)
    {
        FatalErrorIn("createCellZones")
            << "CellZone: " << zoneName << " exists on the mesh."
            << nl << " Check: "
            << "zoneDict"
            << nl << "Solve this by executing:"
            << nl << "> rm constant/polyMesh/*Zones"
            << nl << "> blockMesh"
            << nl << "> createCellZones"            
            << exit(FatalError);
    }    
    
    
    if (dict.found("writeCellSet"))
    {
        writeCellSet_ = Switch(dict.lookup("writeCellSet"));
    }

    const word option(dict.lookup("option"));
    option_ = option;    


    setZone();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

newCellZone::~newCellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void newCellZone::setZone()
{
    DynamicList<label> cells(0);
    
    if (option_ == "boundBox")
    {    
        vector startPoint = dict_.lookup("startPoint");
        vector endPoint = dict_.lookup("endPoint");
    
        boundedBox bb;
        bb.resetBoundedBox(startPoint, endPoint);
        
        for (label i = 0; i < mesh_.nCells(); i++)
        {
            const vector& cellCentreI = mesh_.cellCentres()[i];        
            
            if (bb.contains(cellCentreI))
            {
                cells.append(i);
            }
        }
    }
    
    if (option_ == "elipse")
    {    
        vector centrePoint = dict_.lookup("centrePoint");
        scalar a = readScalar(dict_.lookup("a"));
        scalar b = readScalar(dict_.lookup("b"));
        scalar c = readScalar(dict_.lookup("c"));
        
        for (label i = 0; i < mesh_.nCells(); i++)
        {
            const vector& cellCentreI = mesh_.cellCentres()[i]; 
            
            scalar elipseEquation = sqr(cellCentreI.x() - centrePoint.x())/sqr(a)
                + sqr(cellCentreI.y() - centrePoint.y())/sqr(b)
                + sqr(cellCentreI.z() - centrePoint.z())/sqr(c);
            
            if (elipseEquation <= 1.0)
            {
                cells.append(i);
            }            
        }
    }
    
    if (option_ == "outsideElipse")
    {    
        vector centrePoint = dict_.lookup("centrePoint");
        scalar a = readScalar(dict_.lookup("a"));
        scalar b = readScalar(dict_.lookup("b"));
        scalar c = readScalar(dict_.lookup("c"));
        
        for (label i = 0; i < mesh_.nCells(); i++)
        {
            const vector& cellCentreI = mesh_.cellCentres()[i]; 
            
            scalar elipseEquation = sqr(cellCentreI.x() - centrePoint.x())/sqr(a)
                + sqr(cellCentreI.y() - centrePoint.y())/sqr(b)
                + sqr(cellCentreI.z() - centrePoint.z())/sqr(c);
            
            if (elipseEquation > 1.0)
            {
                cells.append(i);
            }            
        }
    }
    
    cells.shrink();

    cells_.setSize(cells.size());

    forAll(cells, i)
    {
        cells_[i] = cells[i];
    }    

    scalar minX = GREAT;
    scalar maxX = -GREAT;
    scalar minY = GREAT;
    scalar maxY = -GREAT;
    scalar minZ = GREAT;
    scalar maxZ = -GREAT;

    forAll(cells_, c)
    {
        const label& cellI = cells_[c];
        const labelList& pointList = mesh_.cellPoints()[cellI];

        forAll(pointList, p)
        {
            const label& pointI = pointList[p];
            const vector& vPointI = mesh_.points()[pointI];

            if (vPointI.x() < minX)
            {
                minX = vPointI.x();
            }

            if (vPointI.y() < minY)
            {
                minY = vPointI.y();
            }

            if (vPointI.z() < minZ)
            {
                minZ = vPointI.z();
            }

            if (vPointI.x() > maxX)
            {
                maxX = vPointI.x();
            }

            if (vPointI.y() > maxY)
            {
                maxY = vPointI.y();
            }

            if (vPointI.z() > maxZ)
            {
                maxZ = vPointI.z();
            }
        }

    }

    Info << "Dimensions of square: " 
         << mag(maxX - minX) << " x " 
         << mag(maxY - minY) << " x " 
         << mag(maxZ - minZ)
         << endl;
         
}


} // End namespace Foam

// ************************************************************************* //
