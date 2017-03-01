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
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

// --- Default libraries
#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "multivariateScheme.H"
#include "fvIOoptionList.H"

// --- Smoothing for Local-Time-Stepping
#include "fvcSmooth.H"

// --- High-temperature chemistry library
#include "rho2HTCModel.H"

// --- Relaxation processes
#include "relaxationTimeModel.H"
#include "relaxationTimeModelVV.H"
#include "relaxationTimeModelHE.H"
#include "relaxationTimeModeleV.H"

// --- Multi-species transport models
#include "mixingRule.H"
#include "multiSpeciesTransportModel.H"

// --- Rarefied gas library
#include "rarefactionParameter.H"

// --- Numerics
#include "numerics/AUSM_functions.H"

// --- Others
#include "wallFvPatch.H"
#include "constants.H"

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
    
    Info<< "\e[1;33mTotal no of Iterations " << totNoIteration << "\e[0m\n"
        << "\e[1;31mEnd\e[0m\n" << endl;
    
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
    #include "createFields.H"   
    #include "readTimeControls.H"
    
    // --- Initialisation of reacting and two-temperature related fields 
    #include "createReactingFields.H"
    #include "createVibrationalFields.H"
    
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "numerics/readFluxScheme.H"
    
    // --- Write quantities in the 0 folder
    #include "write/write0.H" 

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    /*dimensionedScalar vv_zero("vv_zero", dimArea, 0.0);
    dimensionedVector Vv_one("vv_one", dimArea, vector::one);
    dimensionedScalar vv_one("vv_one", dimArea, 1.0);
    dimensionedScalar vvv_zero("vvv_zero", dimMass/dimTime, 0.0);
    dimensionedScalar vvvv_zero("vvvv_zero", dimless, 0.0);*/
    
    // --- Upwind interpolation of primitive fields on faces
    #include "numerics/upwindInterpolation.H"
    //if(fluxScheme == "AUSM+-up")
    //{
        //#include "numerics/upwindInterpolation_AUSM2.H"
    //}
    //else
    //{
          #include "numerics/upwindInterpolation_KNP.H"
    //}
    
    // --- Time control and local time stepping
    #include "numerics/compressibleCourantNo.H"
    #include "LTS/setInitialrDeltaT.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;
    
    label noSubRestart = 0;
    
    while(runTime.run())
    {
        // --- Re-read run-time modifiable dictionaries
        #include "runTimeEditing/hTCPropertiesDictModification.H"
        #include "runTimeEditing/transportDictModification.H"
        #include "runTimeEditing/twoTemperatureDictModification.H"
        
        // --- Upwind interpolation of primitive fields on faces
        #include "numerics/upwindInterpolation.H"
        //if(fluxScheme == "AUSM+-up")
        //{
            //#include "numerics/upwindInterpolation_AUSM2.H"
        //}
        //else
        //{
            #include "numerics/upwindInterpolation_KNP.H"
        //}
        
        // --- Time control and local time stepping
        #include "numerics/compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "LTS/readLTSTimeControls.H"
        #include "setDeltaT.H"
        if(adjustTimeStep)
        {
            //scalar dtChem = chemistry.solve(runTime.deltaT().value());
            //Info<< "dtChem = " << dtChem << endl;
            //runTime.setDeltaT(min(dtChem, runTime.deltaT().value()));
            //Info<< "Dt = " << runTime.deltaT().value() << endl;
        }

        runTime++;
        Info<< "\e[1;33mTime = " << runTime.timeName() << "\e[0m"<< nl << endl;

        // --- Fields needed for solving equations
        //if(fluxScheme == "AUSM+-up")
        //{
              //#include "numerics/declarations_AUSM2.H"
        //}
        //else
        //{
              #include "numerics/declarations_KNP.H"
        //}

        // --- Local time stepping (LTS)
        if(activateLTS)
        {
            #include "LTS/setrDeltaT.H"
        }
        
        // --- Re-set the switch/counter that serve as a warning if the 
        //     temperature goes unbounded
        thermo.resetTemperatureBoundingInfo();
        
        // --- Solve continuity equation
        #include "eqns/rhoEqn.H"

        // --- Solve species transport and reaction equations
        scalar time = runTime.elapsedCpuTime();
        #include "eqns/YEqn.H"
        Info<< "\e[1;34mYEqn time" << tab << runTime.elapsedCpuTime() - time 
            << " s\e[0m" << endl; 
        
        // --- Solve momentum equations
        #include "eqns/UEqn.H"
        
        if(downgradeSingleT)
        {
            // --- Solve the total energy equation
            //- inviscid
            #include "eqns/eEqnInviscid.H"
            //- viscous
            #include "eqns/eEqnViscous.H"
        }
        else if(downgradeSingleTv)
        {
            // --- Solve the vibrational energy equation
            // --- Solve the total energy equation
            //- inviscid
            #include "eqns/evEqnInviscid_SingleTv.H" 
            #include "eqns/eEqnInviscid.H" 
            //- viscous
            #include "eqns/evEqnViscous_SingleTv.H" 
            #include "eqns/eEqnViscous.H" 
        }
        else if(downgradeSingleVibMode)
        {
            // --- Solve the vibrational energy equations
            // --- Solve the total energy equation
            //- inviscid
            #include "eqns/evEqnInviscid.H" 
            #include "eqns/eEqnInviscid.H"
            //- viscous
            #include "eqns/evEqnViscous.H"
            #include "eqns/eEqnViscous.H"
        }
        else
        {
            // --- Solve the vibrational energy equations (one per vib. mode)
            // --- Solve the total energy equation
            //- inviscid
            //#include "Eqns/evEqn_MultiModes.H"
            //#include "Eqns/eEqnInviscid.H" 
            //- viscous
            //#include "Eqns/eEqnViscous.H"   
        }
        
        // --- Pressure field calculation
        #include "eqns/pEqn.H"
        
        rarefactionParameters().correct(U);
        
        turbulence->correct();
        
        // --- Print a report in the log file if temperature had to be bounded
        thermo.temperatureBoundingReport();

        if(runTime.outputTime())
        {
            runTime.write();
            #include "write/write.H"
        }

        previousIterationTime = 
            max(runTime.elapsedCpuTime()-currentIterationTime, 1e-3);
        Info<< "\e[1;33mPhase no " << noRestart << "." << noSubRestart
            << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << "  Iteration no " << noIteration<<" (" << previousIterationTime 
            << " s)\e[0m" << nl << endl;
        currentIterationTime = runTime.elapsedCpuTime();
        noIteration += 1;
    }

    return false;
}

// ************************************************************************* //
