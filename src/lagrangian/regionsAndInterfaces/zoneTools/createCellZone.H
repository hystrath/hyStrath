void createCellZone
(
    polyMesh& mesh,
    const labelList& cellList,
    const word& name
)
{
    if(cellList.size() > 0)
    {
        if (mesh.cellZones().findZoneID(name) == -1)
        {
            label scz = mesh.cellZones().size();
            mesh.cellZones().setSize(scz+1);
            mesh.cellZones().set
            (
                scz,
                new cellZone
                (
                    name,             //name
                    cellList,                 //addressing
                    scz,                        //index
                    mesh.cellZones()            //cellZoneMesh
                )
            );

            mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
        }
    }
}
