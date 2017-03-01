
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

    scalar timee = runTime.elapsedCpuTime();
    if(speciesDiffusionName != "noSpeciesDiffusion")
    {
        speciesDiffusion().correct();
    }
    Info<< "\e[1;32mcorrectSpeciesDiffusion time" << tab << (runTime.elapsedCpuTime() - timee) << " s\e[0m" << endl;

    reaction->correct();
    dQ = reaction->dQ();

    volScalarField Ytot(0.0*Y[0]);

    for(int speciei=0; speciei < Y.size(); speciei++) // TODO Y.size()-1
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
                - fvc::div(speciesDiffusion().sumDiffusiveFluxes()*Yi, "div(diffusiveFluxes)");
        }
        
        fvOptions.constrain(YiEqn);
        YiEqn.solve(mesh.solver("Yi"));
        fvOptions.correct(Yi);
        
        if(hyDev)
        {
            //speciesDiffusion().getSpeciesMassFlux(speciei, YiEqn.flux()); TODO
        }
        
        Yi.max(0.0);
        Ytot += Yi;
        
        Yi.correctBoundaryConditions();
    }
    
    //Y.last() = min(max(1.0 - Ytot, 0.0), 1.0); // TODO reintroduce
    //Y.last().correctBoundaryConditions();
    // calculate the species mass flux for the last species
}

thermo.correctFractions();

if(hyLight) thermo.correct2T_Light();
else thermo.correct2T();

scalar timee = runTime.elapsedCpuTime();
transport().correct();
Info<< "\e[1;32mcorrectTransport time" << tab << (runTime.elapsedCpuTime() - timee) << " s\e[0m" << endl;

