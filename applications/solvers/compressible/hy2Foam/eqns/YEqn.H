
if(Y.size() > 1 and (speciesDiffusionName != "noSpeciesDiffusion" or reaction->active() or injectors))
{
    tmp<fv::convectionScheme<scalar> > mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );

    if(speciesDiffusionName != "noSpeciesDiffusion")
    {
        scalar timee = runTime.elapsedCpuTime();
        speciesDiffusion().correct();
        Info<< "\e[1;32mcorrectSpeciesDiffusion" << tab << (runTime.elapsedCpuTime() - timee) << " s\e[0m" << endl;
    }
    
    scalar timef = runTime.elapsedCpuTime();
    reaction->correct();
    dQ = reaction->dQ();
    Info<< "\e[1;32mcorrectReactions" << tab << (runTime.elapsedCpuTime() - timef) << " s\e[0m" << endl;

    volScalarField Ytot(0.0*Y[0]);
    
    const volScalarField& Wmix = composition.Wmix();
    
    //for(int speciei=0; speciei < Y.size()-1; speciei++) 
    // TODO two different versions w Y.size()-1
    forAll(Y, speciei) 
    {
        volScalarField& Yi = Y[speciei];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(turbulence->phi(), Yi)
        );
        
        if(reaction->active())
        {
            YiEqn -= reaction->R(Yi);
        }
        
        if(speciesDiffusionName != "noSpeciesDiffusion")
        {
            YiEqn += fvm::laplacian(-speciesDiffusion().rhoD(speciei), Yi, "laplacian(rhoD,Yi)")
                + fvc::div
                  (
                      speciesDiffusion().JGradp(speciei)
                    + speciesDiffusion().JGradT(speciei),
                      "div(diffusiveFluxes)"
                  );
                  
            if(not speciesDiffusion().useNonCorrected())
            {
                YiEqn -= fvc::div(speciesDiffusion().sumDiffusiveFluxes()*Yi, "div(diffusiveFluxes)");
            }      
                
            if(speciesDiffusion().solvingForX())
            {
                YiEqn -= fvc::div(speciesDiffusion().rhoD(speciei)
                    *Yi/Wmix*fvc::grad(Wmix), "div(diffusiveFluxes)");
            }   
        }
        
        fvOptions.constrain(YiEqn);
        YiEqn.solve(mesh.solver("Yi"));
        fvOptions.correct(Yi);
        
        if(hyDev)
        {
            //speciesDiffusion().getSpeciesMassFlux(speciei, YiEqn.flux());
        }
        
        Yi.max(0.0);
        Ytot += Yi;
    }
    
    // TODO reintroduce in the Ns-1 version
    //Y.last() = min(max(1.0 - Ytot, 0.0), 1.0); 
    //Y.last().correctBoundaryConditions();
    // calculate the species mass flux for the last species

    // NB: ensures that SUM(Y) = 1 without risking SUM(Y-1) > 1
    forAll(Y, speciei)
    {
        Y[speciei].internalField() /= Ytot.internalField(); 
    }

    forAll(Y, speciei)
    {
        Y[speciei].correctBoundaryConditions();
    }
}

thermo.correctFractions();

if(hyLight) thermo.correct2T_Light();
else thermo.correct2T();

scalar timee = runTime.elapsedCpuTime();
transport().correct();
Info<< "\e[1;32mcorrectTransport" << tab << (runTime.elapsedCpuTime() - timee) << " s\e[0m" << endl;

