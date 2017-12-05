void cellsToFaceZone
(
    polyMesh& mesh,
    const labelList& cellList,
    const word& name
)
{
    DynamicList<label> faceList(0);

    forAll(cellList, c)
    {
        const label& cellI = cellList[c];

        const labelList& faces = mesh.cells()[cellI];

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
