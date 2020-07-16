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

void createFaceSet
(
    const polyMesh& mesh,
    const labelList& set,
    const word& setName
)
{
//     if(set.size() > 0)
//     {
        // Load set to work
        autoPtr<topoSet> currentSetPtr(NULL);

        //- make a back up of the set
//         backup(mesh, setName, setName+"_old");

        currentSetPtr.reset
        (
            new faceSet
            (
                mesh,
                setName,
                set.size()
            )
        );

        topoSet& currentSet = currentSetPtr();

        // - filling set
        forAll(set, face)
        {
            const label& faceI = set[face];

            currentSet.insert(faceI);
        }

        currentSet.sync(mesh); // removed this cause it was messing about with parallisation

        Info<< "Writing " << currentSet.name()
            << " (size " << currentSet.size() << ") to "
            << currentSet.instance()/currentSet.local()
               /currentSet.name()
            << endl << endl;

        currentSet.write();
//     }
}


// try to read in sets by the following code

//         autoPtr<topoSet> currentSetPtr(NULL);
//         IOobject::readOption r;
//
//         r = IOobject::MUST_READ;
//
//         currentSetPtr.reset
//         (
//             new faceSet
//             (
//                 mesh,
//                 setName,
//                 r
//             )
//         );
