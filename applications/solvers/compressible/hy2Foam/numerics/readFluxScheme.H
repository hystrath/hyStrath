word fluxScheme("Kurganov");
if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if ((fluxScheme == "Tadmor") or (fluxScheme == "Kurganov") or (fluxScheme == "AUSM+-up"))
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else
    {
        FatalErrorIn
        (
            "rhoCentralFoam::readFluxScheme"
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov, AUSM+-up"
            << abort(FatalError);
    }
}

scalar sqrMachInf = sqr(mesh.schemesDict().lookupOrDefault<scalar>("Mach_inf", SMALL));
