//---------------------------------------------------------
// Two-temperature model
//---------------------------------------------------------

// Definition of the relaxation models

const label noSolvedVibEq = composition.solvedVibEqSpecies().size();

// Vibrational-Translational (V-T) energy transfer
autoPtr<relaxationTimeModel> thermalRelaxModelVT;

if (downgradeSingleT)
{
    thermalRelaxModelVT.clear();
}
else
{
    thermalRelaxModelVT = relaxationTimeModel::New
    (
        thermo,
        turbulence()
    );
}

// Vibrational-Vibrational (V-V) energy transfer
autoPtr<relaxationTimeModelVV> thermalRelaxModelVV;

if (downgradeSingleT or downgradeSingleTv or noSolvedVibEq < 2)
{
    thermalRelaxModelVV.clear();
}
else
{
    thermalRelaxModelVV = relaxationTimeModelVV::New
    (
        thermo,
        turbulence()
    );
}

// Heavy-particles - free electrons (h-e) energy transfer
autoPtr<relaxationTimeModelHE> thermalRelaxModelHE;

if (downgradeSingleT or not containsElectrons)
{
    thermalRelaxModelHE.clear();
}
else
{
    thermalRelaxModelHE = relaxationTimeModelHE::New
    (
        thermo,
        turbulence()
    );
}

// Free-electrons - vibrational energy transfer
autoPtr<relaxationTimeModeleV> thermalRelaxModeleV;

if (downgradeSingleT or downgradeSingleTv or not containsElectrons)
{
    thermalRelaxModeleV.clear();
}
else
{
    thermalRelaxModeleV = relaxationTimeModeleV::New
    (
        thermo,
        turbulence()
    );
}


//---------------------------------------------------------
// Loading thermophysical quantities relative to vib. mode
//---------------------------------------------------------

volScalarField rhoEv
(
    IOobject
    (
        "rhoEv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*ev
);

PtrList<volScalarField> rhoEvk(noSolvedVibEq);

if (downgradeSingleT or downgradeSingleTv)
{
   rhoEvk.clear();
}
else
{
    forAll(rhoEvk, speciei)
    {
        rhoEvk.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "rhoEv_" + Y[speciei].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pD[speciei]*hevk[speciei]
            )
        );
    }

    rhoEv *= 0.0;

    forAll(hevk, speciei)
    {
        rhoEv += pD[speciei]*hevk[speciei];
    }

    ev = rhoEv/rho;
}


// ABORTIVE WORK ( multiple vibrational temperatures for triatomic molecules)
/*PtrList<PtrList<volScalarField> >& hevk_mode = composition.hevel_mode();

PtrList<PtrList<volScalarField> > rhoEvk_mode(Y.size());
forAll (rhoEvk_mode, speciei)
{
    rhoEvk_mode.set
    (
        speciei,
        new PtrList<volScalarField>(composition.noVibrationalTemp(speciei))
    );
}

forAll (rhoEvk_mode, speciei)
{
    forAll (rhoEvk_mode[speciei], vibMode)
    {
        rhoEvk_mode[speciei].set
        (
            vibMode,
            new volScalarField
            (
                IOobject
                (
                    "rhoEv_" + Y[speciei].name() + "." + name(vibMode+1),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pD[speciei]*hevk_mode[speciei][vibMode]
            )
        );
    }
}*/
