
// --- Solve vibrational energy equations, VISCOUS
if (!inviscid)
{
    fvScalarMatrix evEqnViscous
    (
        fvm::ddt(rho, ev) - fvc::ddt(rho, ev)
      - fvm::laplacian(transport().alphave(), ev)
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        volVectorField sumJsByevs = speciesDiffusion().Jcorr(0)*hevk[0];

        for(label speciei = 1; speciei < Y.size(); speciei++)
        {
            sumJsByevs += speciesDiffusion().Jcorr(speciei)*hevk[speciei];
        }

        evEqnViscous += fvc::div(sumJsByevs, "div(diffusiveFluxes)");
    }

    evEqnViscous.solve();

    rhoEv = rho*ev;

    thermo.correctHEVK();
}

