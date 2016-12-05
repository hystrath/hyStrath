
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
        speciesDiffusion().correct();
    }

    reaction->correct();
    dQ = reaction->dQ();

    volScalarField Ytot(0.0*Y[0]);

    for(int speciei=0; speciei < Y.size()-1; speciei++)
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
        
        Yi.max(0.0);
        Ytot += Yi;
        
        Yi.correctBoundaryConditions();
    }
    
    Y.last() = min(max(1.0 - Ytot, 0.0), 1.0);
    Y.last().correctBoundaryConditions();
}

thermo.correctFractions();

if(hyLight) thermo.correct2T_Light();
else thermo.correct2T();
#include "eqns/pEqn.H"
transport().correct();

