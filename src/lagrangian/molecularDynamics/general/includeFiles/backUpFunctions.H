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

// Copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const topoSet& fromSet,
    const word& toName
)
{
//     Info<< "Backing up " << fromName << " into " << toName << endl;

    topoSet backupSet(mesh, toName, fromSet);

    backupSet.write();
}


// Read and copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const word& toName
)
{
    topoSet fromSet(mesh, fromName, IOobject::READ_IF_PRESENT);

    backup(mesh, fromName, fromSet, toName);
}
