void pointsToFaceZone
(
    polyMesh& mesh,
    const labelList& pointList,
    const word& name
)
{
    DynamicList<label> faceList(0);

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
        //faceList.shrink(),
		faceList,
        name
    );
}
