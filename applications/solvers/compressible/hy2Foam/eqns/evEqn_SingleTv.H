// The following might be useless VINCENT
thermo.correct();
thermo.correct2T();
transport().correct();

thermalRelaxModelVT().correct();
thermalRelaxModelHE().correct();

rhoEv = rhoEvk[0]*0.0;
surfaceScalarField phiEv = phiEvk[0]*0.0;
volScalarField QVT = thermalRelaxModelVT().QVT(0)*0.0;
volScalarField QCV = reaction->Scv(0)*0.0;

// --- Solve the vibrational energy equation, INVISCID
forAll (rhoEvk, speciei)
{
    rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    rhoEv += rhoEvk[speciei];
    phiEv += phiEvk[speciei];
    
    QVT += thermalRelaxModelVT().QVT(speciei);
    QCV += reaction->Scv(speciei);
}    

solve
(
    fvm::ddt(rhoEv) 
  + fvc::div(phiEv)      
 ==
    QVT
  + thermalRelaxModelHE().QHE()
  + QCV
);

ev = rhoEv/rho;
ev.correctBoundaryConditions();

thermo.correct();
thermo.correct2T();
transport().correct();

rhoEv.boundaryField() = rho.boundaryField()*ev.boundaryField();


// --- Solve the vibrational energy equation, VISCOUS
if (!inviscid)
{
    solve
    (
        fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
      - fvm::laplacian(transport().alphave(), ev)
    );
        
}

rhoEv = rho*ev;

// The following might be useless VINCENT
thermo.correct();
thermo.correct2T();
transport().correct();

