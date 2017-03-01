
// --- Solve vibrational energy equations, VISCOUS
if (!inviscid)
{
    fvScalarMatrix evEqnViscous
    (
        fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
      - fvm::laplacian(transport().alphave(), ev)
    );
    
    if(speciesDiffusionName != "noSpeciesDiffusion")
    {
        forAll(Y, speciei)
        {
            evEqnViscous += 
                fvc::div
                (
                    speciesDiffusion().Jcorrected(speciei)*hevk[speciei],
                    "div(diffusiveFluxes)"
                );
        }
    }
    
    evEqnViscous.solve();
    
    rhoEv = rho*ev;
    
    thermo.correctHEVK();
}

