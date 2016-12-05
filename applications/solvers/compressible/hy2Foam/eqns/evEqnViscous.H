
// --- Solve vibrational energy equations, VISCOUS
if(!inviscid)
{
    forAll(rhoEvk, speciei) // NEW VINCENT 05/08/2016
    {
        fvScalarMatrix evEqnViscous
        (
            fvm::ddt(pD[speciei], hevk[speciei]) - fvc::ddt(pD[speciei], hevk[speciei])
          - fvm::laplacian(transport().alphave(speciei), hevk[speciei])
        );
        
        if(speciesDiffusionName != "noSpeciesDiffusion")
        {
            evEqnViscous += fvc::div(speciesDiffusion().Jcorrected(speciei)*hevk[speciei], "div(diffusiveFluxes)");
        }
        
        evEqnViscous.solve();
        
        rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    }
    
    thermo.correctNonMoleculeHEVK(); // NEW VINCENT 11/08/2016
 
    rhoEv *= 0.0;
    forAll(hevk, speciei) // NEW VINCENT 05/08/2016
    {
        rhoEv += pD[speciei]*hevk[speciei];
    }

    ev = rhoEv/rho;
}

