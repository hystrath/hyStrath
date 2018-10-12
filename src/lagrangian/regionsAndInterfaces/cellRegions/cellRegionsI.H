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
    ovRegion

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// inline const labelList& layerRegions::remainingCells() const
// {
//     return remainingCells_;
// }
// 
// inline const labelList& layerRegions::regionsCells() const
// {
//     return regionsCells_;
// }




/*
inline void ovRegion::setOverlapLayers()
{
    setInitialPointInfo();

    // loop through subRegions and build layers
//     forAll(this*, sR)
//     {
//         this*[sR].setSubLayerPoints(); // build layers
//     }
}*/











/*
inline void ovRegion::setIntermediatePointInfo
(
    const List<pointInfo>& pointData,
    const label& newSize
)
{
    pointDataNew_.clear();
    pointDataNew_.setSize(newSize);
//    List<pointInfo> pointDataTruncated(newSize);

    label count = 0;

    forAll(pointData, p)
    {
	// make a copy of those point infos that belong to the right zone but have not been accepted
	// this should give all those points in the remaining region, including the new interface points

	if(!pointData[p].acceptedPoint())
	{
           pointDataNew_[count] = pointInfo(pointData[p]);
	}
    }
//    return pointDataNew;
}*/



/*
inline void ovRegion::setAcceptedPoints()
{
    DynamicList<label> acceptedPoints(0);

    forAll(pointData_, p)
    {
        if(pointData_[p].acceptedPoint())
        {
            acceptedPoints.append(pointData_[p].currentPointLabel());
        }
    }

    acceptedPoints_ = acceptedPoints.shrink();
}


inline void ovRegion::setNewInterfacePoints()
{
    DynamicList<label> newInterfacePoints(0);

    forAll(pointData_, p)
    {
        if(pointData_[p].pointOnInterface())
        {
            newInterfacePoints.append(pointData_[p].currentPointLabel());
        }
    }

    newInterfacePoints_ = newInterfacePoints.shrink();
}
*/



// - Access

// inline const labelList& ovRegion::acceptedPoints() const
// {
//     return acceptedPoints_;
// }
// 
// 
// inline const labelList& ovRegion::newInterfacePoints() const
// {
//     return newInterfacePoints_;
// }




// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
