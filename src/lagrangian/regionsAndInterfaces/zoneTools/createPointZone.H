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
