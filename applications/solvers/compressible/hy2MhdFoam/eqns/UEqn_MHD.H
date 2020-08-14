
// --- Solve momentum equation, INVISCID
solve(fvm::ddt(rhoU) + fvc::div(phiUp));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
if (hyDev)
{
    thermo.limitVelocityAtWallBoundary(U.boundaryFieldRef());
}

rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();


// --- Solve momentum equation, VISCOUS
if (!inviscid)
{
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(muEff, U)
      - fvc::div(tauMC)
      - mhd->F(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();
    
    fvOptions.constrain(UEqn);
    
    UEqn.solve();
    
    fvOptions.correct(U);
    
    rhoU = rho*U;
}
