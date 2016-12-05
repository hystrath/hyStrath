
// --- Solve energy, INVISCID
surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
      + (mesh.Sf() & fvc::interpolate(tauMC))
    )
    & (a_pos*U_pos + a_neg*U_neg)
);


solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  - fvc::div(sigmaDotU)
 ==
    reaction->Sh()
); 

e = rhoE/rho - 0.5*magSqr(U);

et.dimensionedInternalField() = e.dimensionedInternalField() - ev.dimensionedInternalField();
et.correctBoundaryConditions();
e.boundaryField() = et.boundaryField() + ev.boundaryField();

rhoE.boundaryField() = rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

