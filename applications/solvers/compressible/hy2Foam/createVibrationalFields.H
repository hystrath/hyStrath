//---------------------------------------------------------
// Two-temperature model
//---------------------------------------------------------

// Definition of the relaxation models

label noSolvedVibEq = composition.solvedVibEqSpecies().size(); // NEW VINCENT 05/08/2016

// Vibrational-Translational (V-T) energy transfer
autoPtr<relaxationTimeModel> thermalRelaxModelVT;

if(downgradeSingleT)
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

if(downgradeSingleT or downgradeSingleTv or noSolvedVibEq < 2)
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

if(downgradeSingleT or not containsElectrons)
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

if(downgradeSingleT or downgradeSingleTv or not containsElectrons)
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
    
if(downgradeSingleT) // NEW VINCENT 05-16/08/2016
{
   rhoEvk.clear();
   rhoEv.clear();
}
else if(downgradeSingleTv)
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

