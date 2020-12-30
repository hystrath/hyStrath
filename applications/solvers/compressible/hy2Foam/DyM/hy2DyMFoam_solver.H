// --- Re-read run-time editable dictionaries
#include "runTimeEditing/hTCPropertiesDictModification.H"
#include "runTimeEditing/transportDictModification.H"
#include "runTimeEditing/twoTemperatureDictModification.H"

// --- Time control and local time stepping
#include "numerics/compressibleCourantNo.H"
#include "readTimeControls.H"
#include "LTS/readLTSTimeControls.H"
#include "setDeltaT.H"

if (adjustTimeStep)
{
    //scalar dtChem = chemistry.solve(runTime.deltaT().value());
    //Info<< "dtChem = " << dtChem << endl;
    //runTime.setDeltaT(min(dtChem, runTime.deltaT().value()));
    //Info<< "Dt = " << runTime.deltaT().value() << endl;
}

runTime++;
Info<< "Time = " << runTime.timeName() << nl << endl;

// --- Do any mesh changes
scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

mesh.update();

if (mesh.changing())
{
    noRefinement++;

    Info<< "Execution time for mesh.update() = "
        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
        << " s" << endl;
}

// --- Directed interpolation of primitive fields onto faces
#include "numerics/upwindInterpolation.H"
// --- Compute numerical fluxes
#include "numerics/fluxesCalculation.H"

        
volScalarField muEff(turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

// --- Local time stepping (LTS)
if (activateLTS)
{
    #include "LTS/setrDeltaT.H"
}

// --- Re-set the switch/counter that serves as a warning if the
//     temperature goes unbounded
thermo.resetTemperatureBoundingInfo();

// --- Solve continuity equation
#include "eqns/rhoEqn.H"

// --- Solve species transport and reaction equations
#include "eqns/YEqn.H"

// --- Solve momentum equations
#include HY2FOAM_EXTERNAL_FILE_MOMENTUM_EQNS

if (downgradeSingleT)
{
    // --- Solve the total energy equation
    //- inviscid
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else if (downgradeSingleTv)
{
    // --- Solve the vibrational energy equation
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include "eqns/evEqnViscous_SingleTv.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else if (downgradeSingleVibMode)
{
    // --- Solve the vibrational energy equations
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/evEqnInviscid.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    #include "eqns/evEqnViscous.H"
    #include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}
else
{
    // --- Solve the vibrational energy equations (one per vib. mode)
    // --- Solve the total energy equation
    //- inviscid
    //#include "eqns/evEqn_MultiModes.H"
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYIVC_EQN
    //- viscous
    //#include HY2FOAM_EXTERNAL_FILE_TOTENERGYVIS_EQN
}

// --- Pressure field calculation
#include "eqns/pEqn.H"

// --- Enable input from / coupling with another application
#ifdef HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#include HY2FOAM_EXTERNAL_FILE_HYBRID_COUPLING
#endif

turbulence->correct();

rarefactionParameters().correct(U);


// --- Print a report in the log file if temperature had to be bounded
thermo.temperatureBoundingReport();

// --- Update field to base refinement on
if ( (noIteration + 1)%refineInterval == 0 or runTime.outputTime())
{
    if (useCustomAdapField)
    {
        if (adapFieldName == "normalisedDensityGradient")
        {
            adapField = mag(fvc::grad(rho))/rho;
        }
        else if (adapFieldName == "normalisedPressureGradient")
        {
            adapField = mag(fvc::grad(p))/p;
        }
        else if (adapFieldName == "MachGradient")
        {
            adapField = mag(fvc::grad(Mach));
        }
    }
}

if (runTime.outputTime())
{
    runTime.write();
    #include "write/write.H"

    #ifdef HY2FOAM_EXTERNAL_FILE_OUTPUT
    #include HY2FOAM_EXTERNAL_FILE_OUTPUT
    #endif
}

previousIterationTime =
    max(runTime.elapsedCpuTime()-currentIterationTime, 1e-3);

Info<< "Phase no " << noRefinement << "/" << noRestart << "." << noSubRestart
    << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << "  Iteration no " << noIteration<<" (" << previousIterationTime
    << " s)" << nl << endl;

currentIterationTime = runTime.elapsedCpuTime();
noIteration += 1;
