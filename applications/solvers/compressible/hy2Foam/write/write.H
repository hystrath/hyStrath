
// --- Output at write time
rarefactionParameters().write();
transport().write();
composition.write();

//- Write diffusion fluxes
if (not hyLight)
{
/*    volVectorField sumDiffusiveFluxes = speciesDiffusion().sumDiffusiveFluxes();*/
/*    sumDiffusiveFluxes.write();*/
/*    forAll(Y, speciei)*/
/*    {*/
/*        volVectorField Jnc = speciesDiffusion().J(speciei);*/
/*        Jnc.write();*/

/*        volVectorField Jc*/
/*        (*/
/*            IOobject*/
/*            (*/
/*                "J_" + composition.species()[speciei],*/
/*                runTime.timeName(),*/
/*                mesh,*/
/*                IOobject::NO_READ,*/
/*                IOobject::NO_WRITE*/
/*            ),*/
/*            speciesDiffusion().Jcorr(speciei)*/
/*        );*/
/*        Jc.write();*/
/*    }*/
}

if (not downgradeSingleT)
{
    if (not hyLight)
    {
        thermo.Tov().write();

        if (noSolvedVibEq > 1)
        {
            thermo.Tv().write();
        }
    }
}

if (downgradeSingleTv)
{
    thermo.Tv().write();
}


if (writeWallHeatFlux)
{
    surfaceScalarField heatFlux_convtr =
        -fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(et);

    surfaceScalarField heatFlux_convve =
        -fvc::interpolate(transport().alphave())*fvc::snGrad(ev);
        
    surfaceScalarField heatFlux_diff =
        speciesDiffusion().getDiffusiveWallHeatFlux();

    surfaceScalarField heatFlux = heatFlux_convtr + heatFlux_convve
        + heatFlux_diff;

    const surfaceScalarField::Boundary& patchHeatFlux = heatFlux.boundaryField();

    Info<< "\nIntegrated wall heat flux [W]" << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            wallHeatFlux.boundaryFieldRef()[patchi] = patchHeatFlux[patchi];

            const scalar int_tot_whf =
                gSum
                (
                    mesh.magSf().boundaryField()[patchi]
                   *patchHeatFlux[patchi]
                );
                
            scalar int_convtr_whf = 0.0;
            scalar int_convve_whf = 0.0;
            scalar int_diff_whf = 0.0;    

            if ( mag(int_tot_whf) > SMALL )
            {
                int_convtr_whf = 
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_convtr.boundaryField()[patchi]
                    )/int_tot_whf*100.0;
                    
                int_convve_whf = 
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_convve.boundaryField()[patchi]
                    )/int_tot_whf*100.0;
                    
                int_diff_whf =
                    gSum
                    (
                        mesh.magSf().boundaryField()[patchi]
                       *heatFlux_diff.boundaryField()[patchi]
                    )/int_tot_whf*100.0;  
            }
            
            Info<< "Patch " << patchi
                << " named " << patchHeatFlux[patchi].patch().name()
                << ": " << int_tot_whf << nl
                << "   convective trans-rot  (%): " << int_convtr_whf << nl
                << "   convective vibro-elec (%): " << int_convve_whf << nl
                << "   diffusive             (%): " << int_diff_whf << nl
                << endl;    
        }
    }
    Info<< endl;

    wallHeatFlux.write();
}


if (reaction->lookupOrDefault<bool>("applyChangesAtWriteTime", false))
{
    return true;
}

if (reaction->lookupOrDefault<bool>("applyChangesAtWriteTimeAndWait", false))
{
    label restartIn = reaction->lookupOrDefault<label>("restartIn", 120);

    Info << "Waiting " << restartIn << " seconds before restart\n\n\n";

    Foam::sleep(restartIn);

    return true;
}
