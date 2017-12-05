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
