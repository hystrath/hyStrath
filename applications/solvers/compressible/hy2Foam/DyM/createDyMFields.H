//---------------------------------------------------------
// Creating fields related to mesh motion/refinement
//---------------------------------------------------------

// Read DyM dictionary
IOdictionary refineDict
(
    IOobject
    (
        "dynamicMeshDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

// List of hardcoded fields to base refinement on
wordList hardcodedAdapFields (3, word::null);
hardcodedAdapFields[0] = "normalisedDensityGradient";
hardcodedAdapFields[1] = "normalisedPressureGradient";
hardcodedAdapFields[2] = "MachGradient";

// Field to base refinement on
bool useCustomAdapField = false;

word adapFieldName =
    refineDict.subDict("dynamicRefineFvMeshCoeffs").lookup("field");
    
label refineInterval =
    refineDict.subDict("dynamicRefineFvMeshCoeffs")
              .lookupOrDefault<label>("refineInterval", GREAT);

volScalarField adapField
(
    IOobject
    (
        adapFieldName,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(adapFieldName, dimensionSet(0, -1, 0, 0, 0), 0.0)
);

forAll(hardcodedAdapFields, f)
{
    if (hardcodedAdapFields[f] == adapFieldName)
    {
        useCustomAdapField = true;
        adapField.writeOpt() = IOobject::AUTO_WRITE;
        break;
    }
}
