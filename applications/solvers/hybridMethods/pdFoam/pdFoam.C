/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    pdFoam

Description
    Particle in Cell solver for 3D, transient, multi-species flows

Version 1.1

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pdCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< nl << "Constructing pdCloud " << endl;
    //- initialise logging variables
    int infoCount = 0;

    //- initialise pd Cloud
    pdCloud pd(runTime, "pd", mesh);

    //- calculate initial fields from particle distribution, setup Leapfrog
    Info << "Calculating initial electromagnetic fields:" << endl;
    pd.setupPIC();

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        infoCount++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        pd.evolve();

        //- write out data to run folder
        runTime.write();

        //- check logging interval
        if(infoCount == pd.infoFrequency())
        {
            //calculate parcel data averages
            pd.info();
            infoCount = 0;

            Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
