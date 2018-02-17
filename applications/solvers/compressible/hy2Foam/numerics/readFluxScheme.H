
typedef void (*FnPtr)
(
    const bool&,
    const bool&,
    const vector&,
    const scalar&,
    const scalar&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const vector&,
    const vector&,
    const vector&,
    const vector&,
    const scalar&,
    const scalar&,
    const PtrList<scalar>&,
    const PtrList<scalar>&,
    const scalar&,
    const scalar&,
    scalar&,
    vector&,
    scalar&,
    vector&,
    PtrList<scalar>&,
    scalar&
);

std::map<word, FnPtr> fluxMap;
fluxMap["Kurganov"] = KurganovFlux;
fluxMap["Tadmor"] = TadmorFlux;
fluxMap["KFVS"] = KFVSFlux;
fluxMap["LxF"] = LxFFlux;
fluxMap["ROE"] = ROEFlux;

word fluxScheme("Kurganov");
FixedList<word, 2> blendedSchemes;

if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
    if (
             (fluxScheme == "Tadmor") 
          or (fluxScheme == "Kurganov")
          or (fluxScheme == "AUSM+-up")
          or (fluxScheme == "KFVS")
          or (fluxScheme == "LxF")
          or (fluxScheme == "ROE")
       )
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else if (fluxScheme == "blended")
    {
        if(mesh.schemesDict().readIfPresent("blendedSchemes", blendedSchemes))
        {
            Info<< "fluxScheme for shock vicinity: " << blendedSchemes[0] << endl;
            Info<< "fluxScheme for near-wall region: " << blendedSchemes[1] << endl;
        }
        else
        {
            FatalErrorIn
            (
                "hy2Foam::readFluxScheme"
            )   << "blendedSchemes: " << fluxScheme
                << " is not a valid choice. "
                << "Options are: Kurganov, Tadmor, ROE, KFVS, LxF and AUSM+-up"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "hy2Foam::readFluxScheme"
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Kurganov, Tadmor, ROE, KFVS, LxF, AUSM+ -up and blended"
            << abort(FatalError);
    }
}

scalar sqrMachInf = sqr(mesh.schemesDict().lookupOrDefault<scalar>("Mach_inf", SMALL));


// --- Declarations ---
surfaceScalarField amaxSf
(
    IOobject
    (
        "amaxSf",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("amaxSf", dimVolume/dimTime, 0.0)
);

surfaceVectorField aByU
(
    IOobject
    (
        "aByU",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("aByU", dimVelocity, vector::zero)
);

surfaceVectorField phiUp
(
    IOobject
    (
        "phiUp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("phiUp", dimMass/dimTime*dimVelocity, vector::zero)
);

surfaceScalarField phiEp
(
    IOobject
    (
        "phiEp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("phiEp", dimEnergy/dimTime, 0.0)
);

PtrList<surfaceScalarField> phiEvk(noSolvedVibEq);

if(downgradeSingleT)
{
    phiEvk.clear();
}
else if(downgradeSingleTv)
{
    phiEvk.resize(1);
    
    phiEvk.set
    (
        0,
        new surfaceScalarField
        (
            IOobject
            (
                "phiEvk",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("phiEvk", dimEnergy/dimTime, 0.0)
        )
    );
}
else
{
    forAll(phiEvk, speciei)
    { 
        phiEvk.set
        (
            speciei,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiEvk" + Y[speciei].name(),
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar
                (
                    "phiEvk" + Y[speciei].name(), 
                    dimEnergy/dimTime, 
                    0.0
                )
            )
        );
    }
}

