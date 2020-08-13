
if
(
    (
        speciesDiffusionName != "noSpeciesDiffusion" 
     or reaction->active()
     or injectors
    )
  and
    Y.size() > 1
)
{
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );

    if (speciesDiffusionName != "noSpeciesDiffusion")
    {
        speciesDiffusion().correct();
    }
    
    reaction->correct();
    dQ = reaction->dQ();

    volScalarField Ytot(0.0*Y[0]);

    const volScalarField& Wmix = composition.Wmix();
    
    forAll(Y, speciei) 
    {
        volScalarField& Yi = Y[speciei];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          ==
            fvOptions(rho, Yi)
        );
        
        if (reaction->active())
        {
            YiEqn -= reaction->R(Yi);
        }

        if (speciesDiffusionName != "noSpeciesDiffusion")
        {
            if (containsElectrons and posElectrons == speciei)
            {
                YiEqn += 
                    fvc::div
                    (
                        speciesDiffusion().J(posElectrons),
                        "div(diffusiveFluxes)"
                    );
            }
            else
            {
                YiEqn += 
                    fvm::laplacian
                    (
                       -speciesDiffusion().rhoD(speciei),
                        Yi,
                        "laplacian(rhoD,Yi)"
                    )
                  + fvc::div
                    (
                        speciesDiffusion().JGradp(speciei)
                      + speciesDiffusion().JGradT(speciei),
                        "div(diffusiveFluxes)"
                    );
                    
                if (not speciesDiffusion().useNonCorrected())
                {
                    YiEqn -=
                        fvc::div
                        (
                            speciesDiffusion().sumDiffusiveFluxes()*Yi,
                            "div(diffusiveFluxes)"
                        );
                }      
                    
                if (speciesDiffusion().solvingForX())
                {
                    YiEqn -=
                        fvc::div
                        (
                            speciesDiffusion().rhoD(speciei)*Yi/Wmix
                          * fvc::grad(Wmix),
                            "div(diffusiveFluxes)"
                        );
                }
            }
        }
        
        YiEqn.relax();
        
        fvOptions.constrain(YiEqn);
        
        YiEqn.solve(mesh.solver("Yi"));
        
        fvOptions.correct(Yi);
        
        Yi.max(0.0);
        Ytot += Yi;
    }
    
    forAll(Y, speciei)
    {
        Y[speciei].primitiveFieldRef() /= Ytot.primitiveField();
        Y[speciei].correctBoundaryConditions();
    }
}

thermo.correctFractions();

if (hyLight)
{
    thermo.correct2T_Light();
}
else
{
    thermo.correct2T();
}

transport().correct();
