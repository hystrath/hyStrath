//---------------------------------------------------------
// Local-Time-Stepping (LTS)
//---------------------------------------------------------
IOdictionary fvSchemes
(
    IOobject
    (
        "fvSchemes",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

bool activateLTS = false;

word ddtSchemes = Foam::word(fvSchemes.subDict("ddtSchemes").lookup("default")); 
if(ddtSchemes == "localEuler")
{
    activateLTS = true;
}

volScalarField rDeltaT
(
    IOobject
    (
        "rDeltaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless/dimTime, 1),
    zeroGradientFvPatchScalarField::typeName
);
