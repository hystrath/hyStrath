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

Applications
    hy2MhdFoam: Two-Temperature Open-Source CFD Solver for Hypersonic
             Weakly-Ionised Reacting Flows, with MHD source terms.
             Can be downgraded to a single-temperature model.

Description
    Density-based compressible flow solver

\*---------------------------------------------------------------------------*/

#define HY2FOAM_EXTERNAL_FILE_MOMENTUM_EQNS "eqns/UEqn_MHD.H"
#define HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN "eqns/eEqnInviscid_MHD.H"
#define HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN "eqns/eEqnViscous_MHD.H"
#define HY2FOAM_EXTERNAL_FILE_OUTPUT "write/write_MHD.H"

#include "hy2MhdFoam_include.H"

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
    #include "hy2MhdFoam_createFields.H"

    Info<< "\nStarting time loop\n" << endl;

    label noSubRestart = 0;

    while(runTime.run())
    {
        #include "hy2Foam_solver.H"
        mhd->update();
    }

    return false;
}

// ************************************************************************* //
