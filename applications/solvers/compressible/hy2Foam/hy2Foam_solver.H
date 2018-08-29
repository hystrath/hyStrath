// --- Re-read run-time editable dictionaries
#include "runTimeEditing/hTCPropertiesDictModification.H"
#include "runTimeEditing/transportDictModification.H"
#include "runTimeEditing/twoTemperatureDictModification.H"

// --- Upwind interpolation of primitive fields on faces
#include "numerics/upwindInterpolation.H"
#include "numerics/fluxesCalculation.H"

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

volScalarField muEff(turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

/*if (hyDev)
{
    Zrot = ZrotInf/(1.0 + sqr(pi)/2.0*sqrt(TZ/thermo.Tt()) + (sqr(pi)/4.0 + pi)*TZ/thermo.Tt());
    Z = 3.0*pi/(4.0*(thermo.zetar()+3.0))*Zrot;
    
    mub = mu*(2.0*thermo.zetar()*Z)/(3.0*(thermo.zetar()+3.0));
    bulkToShearRatio = mub/mu;
    
    tauMC = tauMC + mub*tr(fvc::grad(U))*I;
}*/

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
//scalar time = runTime.elapsedCpuTime();
#include "eqns/YEqn.H"
//Info<< "YEqn time" << tab << runTime.elapsedCpuTime() - time << " s" << endl; 

// --- Solve momentum equations
#include "eqns/UEqn.H"

if (downgradeSingleT)
{
    // --- Solve the total energy equation
    //- inviscid
    #include "eqns/eEqnInviscid.H"
    //- viscous
    #include "eqns/eEqnViscous.H"
}
else if (downgradeSingleTv)
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
else if (downgradeSingleVibMode)
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

// --- Enable input from / coupling with another application
#ifdef HY2FOAM_PATH_TO_EXTERNAL_FILE_1
#include HY2FOAM_PATH_TO_EXTERNAL_FILE_1
#endif

turbulence->correct();

rarefactionParameters().correct(U);


// --- Print a report in the log file if temperature had to be bounded
thermo.temperatureBoundingReport();

if (runTime.outputTime())
{
    #ifdef HY2FOAM_PATH_TO_EXTERNAL_FILE_2
    #include HY2FOAM_PATH_TO_EXTERNAL_FILE_2
    #endif
    
    runTime.write();
    #include "write/write.H"
}

previousIterationTime = 
    max(runTime.elapsedCpuTime()-currentIterationTime, 1e-3);
    
Info<< "Phase no " << noRestart << "." << noSubRestart
    << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << "  Iteration no " << noIteration<<" (" << previousIterationTime 
    << " s)" << nl << endl;
    
currentIterationTime = runTime.elapsedCpuTime();
noIteration += 1;
