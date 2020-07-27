
// --- Solve energy, VISCOUS
if (!inviscid)
{
    fvScalarMatrix eEqnViscous
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
      - (  
            fvm::laplacian(turbulence->alphaEff(), e)
          - fvc::laplacian(turbulence->alphaEff(), ev)
        )
      - (
            fvm::laplacian(transport().alphave(), e)
          - fvc::laplacian(transport().alphave(), et)
        )
      ==
        fvOptions(rho, e)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        eEqnViscous += 
            fvc::div
            (
                speciesDiffusion().multiSpeciesHeatSource(),
                "div(diffusiveFluxes)"
            );
    }
    
    eEqnViscous.relax();

    fvOptions.constrain(eEqnViscous);
    
    eEqnViscous.solve();

    fvOptions.correct(e);

    et.ref() = e() - ev();
    et.correctBoundaryConditions();
    e.boundaryFieldRef() = et.boundaryField() + ev.boundaryField();

    rhoE = rho*(e + 0.5*magSqr(U));
}

if (hyLight)
{
    thermo.correct2T_Light();
}
else
{
    thermo.correct2T();
    thermo.correctOverallTemperature();
}
