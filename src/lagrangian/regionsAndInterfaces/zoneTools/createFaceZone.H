void createFaceZone
(
    polyMesh& mesh,
    const labelList& faceList,
    const word& name
)
{
    if(faceList.size() > 0)
    {
        if (mesh.faceZones().findZoneID(name) == -1)
        {
            label sfz = mesh.faceZones().size();
            mesh.faceZones().setSize(sfz+1);
            mesh.faceZones().set
            (
                sfz,
                new faceZone
                (
                    name,                        //name
                    faceList,                 //addressing
                    boolList(faceList.size(), false), //flipMap
                    sfz,                        //index
                    mesh.faceZones()            //cellZoneMesh
                )
            );

            mesh.faceZones().writeOpt() = IOobject::AUTO_WRITE;
        }
    }
}
