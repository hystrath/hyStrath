
// --- Output at write time
rarefactionParameters().write();
transport().write();
composition.write();

if(not downgradeSingleT)
{
    if(noSolvedVibEq > 1)
    {
        thermo.Tv().write();
    }
    if(not hyLight) thermo.T().write();
}


if(hyDev)
{
    //tauMC.write();
    if(heatFluxVector)
    {
        heatFluxVector->write();
    }
}

if(writeHeatFluxes)
{
    if(convectiveHeatFlux and diffusiveHeatFlux)
    {
        *convectiveHeatFlux = turbulence->kappaEff()*fvc::grad(thermo.Tt());
        forAll(Y, speciei)
        {
            *convectiveHeatFlux += transport().kappave(speciei)*fvc::grad(composition.Tv(speciei));
        }

        *diffusiveHeatFlux = speciesDiffusion().multiSpeciesHeatSource();

        convectiveHeatFlux->write();
        diffusiveHeatFlux->write();
    }
}

if(writeWallHeatFlux) 
{       
    surfaceScalarField heatFlux = -fvc::interpolate(turbulence->alphaEff())*fvc::snGrad(et);
    forAll(Y, speciei)
    {
        heatFlux -= fvc::interpolate(transport().alphave(speciei))*fvc::snGrad(hevk[speciei]);
    }   

    const surfaceScalarField::GeometricBoundaryField& patchHeatFlux = heatFlux.boundaryField();

    Info<< "\nWall heat fluxes [W]" << endl;
    forAll(patchHeatFlux, patchi)
    {
        if (isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            Info<< "Patch " << patchi
                << " named " << patchHeatFlux[patchi].patch().name()
                << ": "
                << gSum
                   (
                       mesh.magSf().boundaryField()[patchi]
                      *patchHeatFlux[patchi]
                   )
                << endl;
        }
    }
    
    Info<< endl;

    volScalarField wallHeatFlux
    (
        IOobject
        (
            "wallHeatFlux",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("wallHeatFlux", dimensionSet(1, 0, -3, 0, 0), 0.0)
    );

    forAll(wallHeatFlux.boundaryField(), patchi)
    {
        wallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
    }

    wallHeatFlux.write();
}


if(reaction->lookupOrDefault<bool>("applyChangesAtWriteTime", false))
{
    return true;
}

if(reaction->lookupOrDefault<bool>("applyChangesAtWriteTimeAndWait", false))
{
    label restartIn = reaction->lookupOrDefault<label>("restartIn", 120);
    
    Info << "Waiting " << restartIn << " seconds before restart\n\n\n";
         
    Foam::sleep(restartIn);     
    
    return true;
}

