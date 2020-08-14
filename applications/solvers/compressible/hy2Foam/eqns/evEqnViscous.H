
// --- Solve vibrational energy equations, VISCOUS
if (!inviscid)
{
    forAll(rhoEvk, speciei)
    {
        DynamicList<label> removedCellsFromMatrix;
        DynamicList<scalar> valuesToImpose;

        forAll(Y[speciei], celli)
        {
            if
            (
                Y[speciei][celli]
                    < Foam::rho2ReactionThermo::miniYforSolvingEvEqn
            )
            {
                removedCellsFromMatrix.append(celli);
            }
        }

        forAll(removedCellsFromMatrix, celli)
        {
            valuesToImpose.append(0.0);
        }

        fvScalarMatrix evEqnViscous
        (
            fvm::ddt(pD[speciei], hevk[speciei])
          - fvc::ddt(pD[speciei], hevk[speciei])
          - fvm::laplacian(transport().alphave(speciei), hevk[speciei])
        );

        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            evEqnViscous +=
                fvc::div
                (
                    speciesDiffusion().Jcorr(speciei)*hevk[speciei],
                    "div(diffusiveFluxes)"
                );
        }

        evEqnViscous.setValues(removedCellsFromMatrix, valuesToImpose);

        evEqnViscous.solve();

        rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    }

    thermo.correctHEVK();

    rhoEv *= 0.0;
    forAll(hevk, speciei)
    {
        rhoEv += pD[speciei]*hevk[speciei];
    }

    ev = rhoEv/rho;
}

