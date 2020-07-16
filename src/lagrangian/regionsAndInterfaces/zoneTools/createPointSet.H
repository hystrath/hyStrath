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

void createPointSet
(
    const polyMesh& mesh,
    const labelList& set,
    const word& setName
)
{
    // Load set to work
    autoPtr<topoSet> currentSetPtr(NULL);

    //- make a back up of the set
//     backup(mesh, setName, setName + "_old");

    currentSetPtr.reset
    (
        new pointSet
        (
            mesh,
            setName,
            set.size()
        )
    );

    topoSet& currentSet = currentSetPtr();

    // - filling set

    forAll(set, point)
    {
        const label& pointI = set[point];

        currentSet.insert(pointI);
    }

    currentSet.sync(mesh);

//     Info<< "Writing " << currentSet.name()
//         << " (size " << currentSet.size() << ") to "
//         << currentSet.instance()/currentSet.local()
//             /currentSet.name()
//         << endl << endl;

    currentSet.write();
}
