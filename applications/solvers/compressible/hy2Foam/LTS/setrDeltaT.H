/*Info<< "Time scales min/max:" << endl;

// Cache old reciprocal time scale field
volScalarField rDeltaT0("rDeltaT0", rDeltaT);

// Flow time scale
{
    rDeltaT.dimensionedInternalField() =
    (
        fvc::surfaceSum(mag(phi))().dimensionedInternalField()
       /((2*maxCo)*mesh.V()*rho.dimensionedInternalField())
    );

    //---------------------------------------------------------
    // Limit the largest time scale
    //---------------------------------------------------------
    rDeltaT.max(1/maxDeltaT);

    Info<< "    Flow        = "
        << gMin(1/rDeltaT.internalField()) << ", "
        << gMax(1/rDeltaT.internalField()) << endl;
}

// Reaction source time scale
if (alphaTemp < 1.0)
{
    volScalarField::DimensionedInternalField rDeltaTT
    (
        mag(reaction->Sh())/(alphaTemp*rho*thermo.Cv_t()*thermo.T()) // NEW VINCENT: was Cp*T
    );


    Info<< "    Temperature = "
        << gMin(1/(rDeltaTT.field() + VSMALL)) << ", "
        << gMax(1/(rDeltaTT.field() + VSMALL)) << endl;

    rDeltaT.dimensionedInternalField() = max
    (
        rDeltaT.dimensionedInternalField(),
        rDeltaTT
    );
}

// Update tho boundary values of the reciprocal time-step
rDeltaT.correctBoundaryConditions();

// Spatially smooth the time scale field
if (rDeltaTSmoothingCoeff < 1.0)
{
    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
}

// Limit rate of change of time scale
// - reduce as much as required
// - only increase at a fraction of old time scale
if
(
    rDeltaTDampingCoeff < 1.0
 && runTime.timeIndex() > runTime.startTimeIndex() + 1
)
{
    rDeltaT = max
    (
        rDeltaT,
        (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
    );
}

Info<< "    Overall     = "
    << gMin(1/rDeltaT.internalField())
    << ", " << gMax(1/rDeltaT.internalField()) << endl;


*/ //TODO VINCENT :: UNSTABLE

///*****************************************************************************************

//---------------------------------------------------------
// Set the reciprocal time-step from the local Courant number
//---------------------------------------------------------
rDeltaT.ref() = max
(
  1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
  fvc::surfaceSum(amaxSf)()()
          /(maxCo*mesh.V())
);

Info<< "    Flow        = "
        << gMin(1/rDeltaT.primitiveField()) << ", "
        << gMax(1/rDeltaT.primitiveField()) << endl;


//---------------------------------------------------------
// Limit the largest time scale
//---------------------------------------------------------
rDeltaT.max(1/maxDeltaT);


//---------------------------------------------------------
// Spatially smooth the time scale field
//---------------------------------------------------------
scalar rDeltaTSmoothingCoeff = 0.1;
fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

Info<< "    Overall     = " << min(1/rDeltaT).value()
    << ", " << max(1/rDeltaT).value() << nl << endl;


