void createPointZone
(
    polyMesh& mesh,
    const labelList& pointList,
    const word& name
)
{
    if (mesh.pointZones().findZoneID(name) == -1)
    {
        label spz = mesh.pointZones().size();
        mesh.pointZones().setSize(spz+1);
        mesh.pointZones().set
        (
            spz,
            new pointZone
            (
                name,             //name
                pointList,                 //addressing
                spz,                        //index
                mesh.pointZones()            //cellZoneMesh
            )
        );

        mesh.pointZones().writeOpt() = IOobject::AUTO_WRITE;
    }
}
