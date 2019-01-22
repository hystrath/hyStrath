
// --- Solve momentum equation, INVISCID
solve(fvm::ddt(rhoU) + fvc::div(phiUp));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
if(hyDev) thermo.limitVelocityAtWallBoundary(U.boundaryFieldRef()); // NEW VINCENT 14/09/2016

rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();


// --- Solve momentum equation, VISCOUS
if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(muEff, U)
      - fvc::div(tauMC)
      - mhd->F(U)
    );
    
    rhoU = rho*U;
}
