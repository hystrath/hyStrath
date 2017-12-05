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
    randomPointsInSquareZone

Description

\*----------------------------------------------------------------------------*/

#include "randomPointsInSquareZone.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void randomPointsInSquareZone:: setInitialConfiguration
(
    const polyMesh& mesh,
    const label& regionId            
)
{
    const cellZoneMesh& cellZones = mesh.cellZones();

    //-set the total volume
    const labelList& cells = cellZones[regionId];

    scalar XMin = GREAT;
    scalar YMin = GREAT;
    scalar ZMin = GREAT;

    scalar XMax = -GREAT;
    scalar YMax = -GREAT;
    scalar ZMax = -GREAT;

    forAll(cells, c)
    {
//         zoneOnMesh = true;

        const label& cellI = cells[c];
        const labelList& cP = mesh.cellPoints()[cellI];

        forAll(cP, p)
        {
            const label& pI = cP[p];
            const vector& pointI = mesh.points()[pI];

            if(pointI.x() < XMin)
            {
                XMin = pointI.x();
            }
            if(pointI.x() > XMax)
            {
                XMax = pointI.x();
            }


            if(pointI.y() < YMin)
            {
                YMin = pointI.y();
            }
            if(pointI.y() > YMax)
            {
                YMax = pointI.y();
            }

            if(pointI.z() < ZMin)
            {
                ZMin = pointI.z();
            }
            if(pointI.z() > ZMax)
            {
                ZMax = pointI.z();
            }
        }
    }

    //- parallel communication
    if(Pstream::parRun())
    {
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << XMin << XMax << YMin
                                 << YMax << ZMin << ZMax /*<< zoneOnMesh*/;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalar XMinP;
                scalar YMinP;
                scalar ZMinP;
            
                scalar XMaxP;
                scalar YMaxP;
                scalar ZMaxP;
//                 bool zoneOnMeshP;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> XMinP >> XMaxP >> YMinP 
                                    >> YMaxP >> ZMinP >> ZMaxP /*>> zoneOnMeshP*/;
                }

//                 if(zoneOnMeshP)
//                 {
                    if(XMinP < XMin)
                    {
                        XMin = XMinP;
                    }
                    if(XMaxP > XMax)
                    {
                        XMax = XMaxP;
                    }
        
                    if(YMinP < YMin)
                    {
                        YMin = YMinP;
                    }
                    if(YMaxP > YMax)
                    {
                        YMax = YMaxP;
                    }
        
                    if(ZMinP < ZMin)
                    {
                        ZMin = ZMinP;
                    }
                    if(ZMaxP > ZMax)
                    {
                        ZMax = ZMaxP;
                    }
//                 }
            }
        }
    }

    XMin_ = XMin;
    YMin_ = YMin;
    ZMin_ = ZMin;

    deltaX_ = XMax - XMin;
    deltaY_ = YMax - YMin;
    deltaZ_ = ZMax - ZMin;

    Info << nl << "zone information: " << endl;

    Info<< " XMin: " << XMin_ << " YMin: " << YMin_ << " ZMin: " << ZMin_
        << " XMax: " << XMax << " YMax: " << YMax << " ZMax: " << ZMax
        << " deltaX: " << deltaX_  << " deltaY: " <<  deltaY_ << " deltaZ: " << deltaZ_
        << endl;

}






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
randomPointsInSquareZone::randomPointsInSquareZone()
:
    deltaX_(0.0),
    deltaY_(0.0),
    deltaZ_(0.0),
    XMin_(0.0),
    YMin_(0.0),
    ZMin_(0.0),
    rndGen_(clock::getTime())
{}



randomPointsInSquareZone::randomPointsInSquareZone
(
    const polyMesh& mesh,
    const label& regionId
)
:
    deltaX_(0.0),
    deltaY_(0.0),
    deltaZ_(0.0),
    XMin_(0.0),
    YMin_(0.0),
    ZMin_(0.0),
    rndGen_(clock::getTime())
{
    Info << "randomPointsInSquareZone: picking points randomly from a square region."
         << endl;

    setInitialConfiguration
    (
        mesh,
        regionId            
    );
}

randomPointsInSquareZone::randomPointsInSquareZone
(
    const boundBox& bb
)
:
    deltaX_(0.0),
    deltaY_(0.0),
    deltaZ_(0.0),
    XMin_(0.0),
    YMin_(0.0),
    ZMin_(0.0),
    rndGen_(clock::getTime())
{
    Info << "randomPointsInSquareZone: picking points randomly from a square region."
         << endl;

    setInitialConfiguration
    (
        bb          
    );
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

randomPointsInSquareZone::~randomPointsInSquareZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void randomPointsInSquareZone::setInitialConfiguration
(
    const boundBox& bb
)
{
    XMin_ = bb.min().x();
    YMin_ = bb.min().y();
    ZMin_ = bb.min().z();

    deltaX_ = bb.span().x();
    deltaY_ = bb.span().y();
    deltaZ_ = bb.span().z();

    Info << nl << "zone information: " << endl;

    Info<< " XMin: " << XMin_ << " YMin: " << YMin_ << " ZMin: " << ZMin_
        << " deltaX: " << deltaX_  << " deltaY: " <<  deltaY_ << " deltaZ: " << deltaZ_
        << endl;
}

vector randomPointsInSquareZone::randomPoint()
{
    vector pointI
    (
        deltaX_*rndGen_.scalar01() + XMin_,
        deltaY_*rndGen_.scalar01() + YMin_,
        deltaZ_*rndGen_.scalar01() + ZMin_
    );

    return pointI;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
