// --- upwind interpolation of primitive fields on faces
surfaceScalarField rho_pos
(
    "rho_pos",
    fvc::interpolate(rho, pos, "reconstruct(rho)")
);
surfaceScalarField rho_neg
(
    "rho_neg",
    fvc::interpolate(rho, neg, "reconstruct(rho)")
);

// NEW VINCENT ********************************************07/11/2015**********
PtrList<surfaceScalarField> pD_pos(noSolvedVibEq);
PtrList<surfaceScalarField> pD_neg(noSolvedVibEq);

if(downgradeSingleT or downgradeSingleTv) // NEW VINCENT 05-16/08/2016
{
    pD_pos.clear();
    pD_neg.clear();
}
else
{
    forAll(pD_pos, speciei)
    { 
        pD_pos.set
        (
            speciei,
            new surfaceScalarField
            (
                "pD_pos_" + Y[speciei].name(),
                fvc::interpolate(pD[speciei], pos, "reconstruct(rho)")
            )
        );
    }

    forAll(pD_neg, speciei)
    { 
        pD_neg.set
        (
            speciei,
            new surfaceScalarField
            (
                "pD_neg_" + Y[speciei].name(),
                fvc::interpolate(pD[speciei], neg, "reconstruct(rho)")
            )
        );
    }
}

// END NEW VINCENT ****************************************07/11/2015**********

surfaceVectorField rhoU_pos
(
    "rhoU_pos",
    fvc::interpolate(rhoU, pos, "reconstruct(U)")
);
surfaceVectorField rhoU_neg
(
    "rhoU_neg",
    fvc::interpolate(rhoU, neg, "reconstruct(U)")
);

volScalarField rPsi(1.0/psi);

surfaceScalarField rPsi_pos
(
    "rPsi_pos",
    fvc::interpolate(rPsi, pos, "reconstruct(T)")
);
surfaceScalarField rPsi_neg
(
    "rPsi_neg",
    fvc::interpolate(rPsi, neg, "reconstruct(T)")
);

surfaceScalarField e_pos
(
    "e_pos",
    fvc::interpolate(e, pos, "reconstruct(T)")
);
surfaceScalarField e_neg
(
    "e_neg",
    fvc::interpolate(e, neg, "reconstruct(T)")
);

// NEW VINCENT ****************************************************************
PtrList<surfaceScalarField> evk_pos(noSolvedVibEq);
PtrList<surfaceScalarField> evk_neg(noSolvedVibEq);

if(downgradeSingleT) // NEW VINCENT 05-16/08/2016
{
    evk_pos.clear();
    evk_neg.clear();  
}
else if(downgradeSingleTv)
{
    evk_pos.resize(1);
    evk_neg.resize(1);
    
    evk_pos.set
    (
        0,
        new surfaceScalarField
        (
            "ev_pos",
            fvc::interpolate(ev, pos, "reconstruct(T)")
        )
    );
    
    evk_neg.set
    (
        0,
        new surfaceScalarField
        (
            "ev_neg",
            fvc::interpolate(ev, neg, "reconstruct(T)")
        )
    );
}
else
{
    forAll(evk_pos, speciei)
    { 
        evk_pos.set
        (
            speciei,
            new surfaceScalarField
            (
                "evk_pos_" + Y[speciei].name(),
                fvc::interpolate(hevk[speciei], pos, "reconstruct(T)")
            )
        );
    }

    forAll(evk_neg, speciei)
    { 
        evk_neg.set
        (
            speciei,
            new surfaceScalarField
            (
                "evk_neg_" + Y[speciei].name(),
                fvc::interpolate(hevk[speciei], neg, "reconstruct(T)")
            )
        );
    }
}
// END NEW VINCENT ************************************************************

volScalarField gammatr("gammatr", thermo.Cp_t()/thermo.Cv_t()); // NEW VINCENT (also useful for bdry cdts)
volScalarField c(sqrt(gammatr*rPsi)); // NEW VINCENT (see Candler 2009)
Mach = mag(U)/c;

surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

