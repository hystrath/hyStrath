/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Applications
    hyFoam: Single-Temperature Open-Source CFD Solver for Supersonic 
            Combusting Flows
             
    hy2Foam: Two-Temperature Open-Source CFD Solver for Hypersonic 
             Weakly-Ionised Reacting Flows

Description
    Density-based compressible flow solver

\*---------------------------------------------------------------------------*/

#include "hy2Foam_include.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool run
(
    argList& args,
    Time& runTime,
    fvMesh& mesh,
    scalar& currentIterationTime,
    scalar& previousIterationTime,
    label& noRestart,
    label& noIteration
);

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    scalar currentIterationTime = 0.0;
    scalar previousIterationTime = 1.0;
    label noRestart = 0;
    label noIteration = 1;
    label totNoIteration = 0;
    bool restart = false;
    
    do
    {
        noRestart += 1;
        
        restart = run
        (
            args,
            runTime,
            mesh,
            currentIterationTime, 
            previousIterationTime,
            noRestart, 
            noIteration
        );
        
        totNoIteration += noIteration - 1;
        noIteration = 1;
        
    } while(restart);
    
    Foam::rho2ReactionThermo::hasCrashedButRecoveredReport();
    
    Info<< "Total no of Iterations " << totNoIteration << "\n"
        << "End\n" << endl;
    
    return 0;
}

bool run
(
    argList& args,
    Time& runTime, 
    fvMesh& mesh,
    scalar& currentIterationTime, 
    scalar& previousIterationTime,
    label& noRestart,
    label& noIteration
)
{
    #include "hy2Foam_createFields.H"

    Info<< "\nStarting time loop\n" << endl;
    
    label noSubRestart = 0;
    
    while(runTime.run())
    {
        #include "hy2Foam_solver.H"
    }

    return false;
}

// ************************************************************************* //
