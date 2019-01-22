    
// --- Solve energy, VISCOUS
if(!inviscid)
{
    fvScalarMatrix eEqnViscous
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
      - fvm::laplacian(turbulence->alphaEff(), e) + fvc::laplacian(turbulence->alphaEff(), ev)
      - fvm::laplacian(transport().alphave(), e) + fvc::laplacian(transport().alphave(), et)
 //     - Q_mag //--------------Alexey
    );
     Info<< "Q_mag is working" <<endl;
    if(speciesDiffusionName != "noSpeciesDiffusion")
    {
        eEqnViscous += fvc::div(speciesDiffusion().multiSpeciesHeatSource(), "div(diffusiveFluxes)");
    }
    
    eEqnViscous.solve();
    
    et.ref() = e() - ev();
    et.correctBoundaryConditions();
    e.boundaryFieldRef() = et.boundaryField() + ev.boundaryField();
    
    rhoE = rho*(e + 0.5*magSqr(U));
}

//scalar timer = runTime.elapsedCpuTime();
if(hyLight) thermo.correct2T_Light();
else
{
    thermo.correct2T();
    thermo.correctOverallTemperature();
}
//Info<< "correctThermo" << tab << (runTime.elapsedCpuTime() - timer) << " s" << endl;

