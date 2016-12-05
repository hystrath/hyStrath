
// --- Solve momentum equation, INVISCID
solve(fvm::ddt(rhoU) + fvc::div(phiUp));

U.dimensionedInternalField() = rhoU.dimensionedInternalField()/rho.dimensionedInternalField();
U.correctBoundaryConditions();
if(hyDev) thermo.limitVelocityAtWallBoundary(U.boundaryField()); // NEW VINCENT 14/09/2016

rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();


// --- Solve momentum equation, VISCOUS
if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(muEff, U)
      - fvc::div(tauMC)
    );
    
    rhoU = rho*U;
}
