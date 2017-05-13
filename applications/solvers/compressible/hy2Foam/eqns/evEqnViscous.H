
// --- Solve vibrational energy equations, VISCOUS
if(!inviscid)
{
    forAll(rhoEvk, speciei)
    {
        // NEW VINCENT 29/12/2016 *********************************************
        DynamicList<label> removedCellsFromMatrix;
        DynamicList<scalar> valuesToImpose;
        
        forAll(Y[speciei], celli)
        {
            if(Y[speciei][celli] < Foam::rho2ReactionThermo::miniYforSolvingEvEqn)
            {
                removedCellsFromMatrix.append(celli);
            }
        }
        
        forAll(removedCellsFromMatrix, celli)
        {
            valuesToImpose.append(0.0);
        }
        // END NEW VINCENT 29/12/2016 *****************************************
        
        fvScalarMatrix evEqnViscous
        (
            fvm::ddt(pD[speciei], hevk[speciei]) - fvc::ddt(pD[speciei], hevk[speciei])
          - fvm::laplacian(transport().alphave(speciei), hevk[speciei])
        );
        
        if(speciesDiffusionName != "noSpeciesDiffusion")
        {
            evEqnViscous += 
                fvc::div
                (
                    speciesDiffusion().Jcorr(speciei)*hevk[speciei],
                    "div(diffusiveFluxes)"
                );
        }
        
        evEqnViscous.setValues(removedCellsFromMatrix, valuesToImpose); // NEW VINCENT 29/12/2016
        
        evEqnViscous.solve();
        
        rhoEvk[speciei] = pD[speciei]*hevk[speciei];
    }
    
    thermo.correctHEVK(); // NEW VINCENT 11/08/2016
 
    rhoEv *= 0.0;
    forAll(hevk, speciei)
    {
        rhoEv += pD[speciei]*hevk[speciei];
    }

    ev = rhoEv/rho;
}

