/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

Applications
    hyDyMFoam: Single-Temperature Open-Source CFD Solver for Supersonic 
            Combusting Flows
             
    hy2DyMFoam: Two-Temperature Open-Source CFD Solver for Hypersonic 
             Weakly-Ionised Reacting Flows
             
Group
    grpCompressibleSolvers grpMovingMeshSolvers             

Description
    Density-based compressible flow solver with support for mesh-motion and 
    topology changes.

\*---------------------------------------------------------------------------*/

#define HY2FOAM_EXTERNAL_FILE_MOMENTUM_EQNS "eqns/UEqn.H"
#define HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN "eqns/eEqnInviscid.H"
#define HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN "eqns/eEqnViscous.H"

#include "hy2Foam_include.H"
#include "dynamicFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool run
(
    argList& args,
    Time& runTime,
    dynamicFvMesh& mesh,
    scalar& currentIterationTime,
    scalar& previousIterationTime,
    label& noRefinement,
    label& noRestart,
    label& noIteration
);

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    scalar currentIterationTime = 0.0;
    scalar previousIterationTime = 1.0;
    label noRefinement = 0;
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
            noRefinement,
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
    dynamicFvMesh& mesh,
    scalar& currentIterationTime,
    scalar& previousIterationTime,
    label& noRefinement,
    label& noRestart,
    label& noIteration
)
{
    #include "hy2Foam_createFields.H"
    #include "createDyMFields.H"

    Info<< "\nStarting time loop\n" << endl;

    label noSubRestart = 0;

    while(runTime.run())
    {
        #include "hy2DyMFoam_solver.H"
    }

    return false;
}

// ************************************************************************* //
