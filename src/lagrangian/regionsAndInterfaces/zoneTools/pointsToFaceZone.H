void pointsToFaceZone
(
    polyMesh& mesh,
    const labelList& pointList,
    const word& name
)
{
    DynamicList<label> faceList(0);
/*
    forAll(pointList, p)
    {
        const labelList& faces = mesh.pointFaces()[pointList[p]];

        forAll(faces, f)
        {
            const label& faceI = faces[f];

            if(findIndex(faceList, faceI) == -1)
            {
                label noOfPoints = 0;

                const labelList& points = mesh.faces()[faceI];

                forAll(points, pS)
                {
                    if(findIndex(pointList, points[pS]) != -1)
                    {
                        noOfPoints++;
                    }
                }

                if(noOfPoints == points.size())
                {
                    faceList.append(faceI);
                }
            }
        }
    }

    faceList.shrink();


    // - check for missed points

    DynamicList<label> acceptedPoints(0);

    forAll(faceList, f)
    {
        const labelList& points = mesh.faces()[faceList[f]];

        forAll(points, p)
        {
            if(findIndex(pointList, points[p]) != -1)
            {
                acceptedPoints.append(points[p]);
            }
        }
    }

    acceptedPoints.shrink();

    DynamicList<label> missedPoints(0);

    forAll(pointList, p)
    {
        if(findIndex(acceptedPoints, pointList[p]) == -1)
        {
            missedPoints.append(pointList[p]);
        }
    }


*/
//     Info << "number of missed points: " << missedPoints.shrink() << endl;


    forAll(pointList, p)
    {
        const labelList& faces = mesh.pointFaces()[pointList[p]];

        forAll(faces, f)
        {
            if(findIndex(faceList, faces[f]) == -1)
            {
                faceList.append(faces[f]);
            }
        }
    }



    createFaceZone
    (
        mesh,
        faceList.shrink(),
        name
    );
}
