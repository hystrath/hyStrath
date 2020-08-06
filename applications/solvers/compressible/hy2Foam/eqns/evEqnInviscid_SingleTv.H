
// --- Solve vibrational energy equation, INVISCID
thermalRelaxModelVT().correct();
volScalarField QVT = thermalRelaxModelVT().QVT();

if (containsElectrons)
{
    thermalRelaxModelHE().correct();
    QVT += thermalRelaxModelHE().QHE();
}

fvScalarMatrix rhoEvInviscid
(
    fvm::ddt(rhoEv)
  + fvc::div(phiEvk[0])
  - QVT
  - reaction->Scv()
  - reaction->Seiir()
);

if (containsElectrons)
{
    rhoEvInviscid += thermo.pe()*fvc::div(U);
}

rhoEvInviscid.solve
(
    mesh.solver("rhoEv")
);

ev.correctBoundaryConditions();
rhoEv.boundaryFieldRef() = rho.boundaryField()*ev.boundaryField();

thermo.correctHEVK();

ev = rhoEv/rho;

