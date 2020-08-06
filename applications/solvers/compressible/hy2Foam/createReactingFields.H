//---------------------------------------------------------
// Creating fields related to reacting mixtures
//---------------------------------------------------------

// Mixture fields
PtrList<volScalarField>& Y = composition.Y();
PtrList<volScalarField>& pD = composition.pD();

Switch containsElectrons = false;
label posElectrons = -1;
if (composition.contains("e-"))
{
    containsElectrons = true;
    posElectrons = composition.species()["e-"];
}

/*Generic multi-variate discretisation scheme class which may be
instantiated for any of the NVD, CNVD or NVDV schemes.*/
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
    fields.add(pD[i]); // NEW VINCENT 07/11/2015
    fields.add(hevk[i]); // NEW VINCENT 10/08/2016
}
fields.add(e);


// Energy release per unit time
volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

//---------------------------------------------------------
// Definition of the multi-species transport model
//---------------------------------------------------------

autoPtr<multiSpeciesTransportModel> speciesDiffusion
(
    multiSpeciesTransportModel::New
    (
        thermo,
        turbulence()
    )
);

word speciesDiffusionName =
    thermo.transportDictionary().subDict("transportModels")
        .lookupOrDefault<word>("multiSpeciesTransport", "noSpeciesDiffusion");

