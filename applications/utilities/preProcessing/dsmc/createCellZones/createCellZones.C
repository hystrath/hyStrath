/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application


Description
    Creates a series of user defined square zones on the mesh.


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "Time.H"
#include "polyMesh.H"


#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

#include "newCellZone.H"

using namespace Foam;

#include "createCellSet.H"
#include "createCellZone.H"
#include "createFaceZone.H"
#include "cellsToFaceZone.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "readZonesDict.H"

    forAll(regionList, sL)
    {

        Info << sL << "  creating cell zone: " << regionZoneNames[sL] << endl;

        newCellZone sZ
        (
            mesh,
            dictionaries[sL]
        );

        if (sZ.cells().size() == 0)
        {
            Info << "WARNING: no cells have been created -- check" << endl;
        }
        
        Info << " number of cells: " << sZ.cells().size() << nl<< endl;
        
        createCellZone(mesh, sZ.cells(), sZ.name());
    
        // visualisation of region
        cellsToFaceZone(mesh, sZ.cells(), sZ.name());

        if (sZ.writeCellSet())
        {
            createCellSet(mesh, sZ.cells(), sZ.name());
        }
        
        Info << endl;
    }

    if (!mesh.write())
    {
            FatalErrorIn(args.executable())
                << "Failed writing cellZones."
                << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl; 

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
