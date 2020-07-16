/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

void pointsToCellsToFaceZone
(
    polyMesh& mesh,
    const labelList& pointList,
    const word& name
)
{
    DynamicList<label> cellList(0);

    forAll(pointList, p)
    {
        const labelList& cells = mesh.pointCells()[pointList[p]];

        forAll(cells, c)
        {
            const label& cellI = cells[c];

            if(findIndex(cellList, cellI) == -1)
            {
                label noOfPoints = 0;

                const labelList& points = mesh.cellPoints()[cellI];

                forAll(points, pS)
                {
                    if(findIndex(pointList, points[pS]) != -1)
                    {
                        noOfPoints++;
                    }
                }

                if(noOfPoints == points.size())
                {
                    cellList.append(cellI);
                }
            }
        }
    }

//     Info << "Accepted Points: "<< pointList << endl;
//
//     Info << "Accepted Cells: "<< cellList.shrink() << endl;

    // check that all points in pointList are included at least once in the cellList.

    cellsToFaceZone
    (
        mesh,
        //cellList.shrink(),
		cellList,
        name
    );
}
