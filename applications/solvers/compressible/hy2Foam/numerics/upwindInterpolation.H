// --- Directed interpolation of primitive fields onto faces

surfaceScalarField rho_pos(interpolate(rho, pos));
surfaceScalarField rho_neg(interpolate(rho, neg));

PtrList<surfaceScalarField> pD_pos(noSolvedVibEq);
PtrList<surfaceScalarField> pD_neg(noSolvedVibEq);

if (downgradeSingleT or downgradeSingleTv)
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
                interpolate(pD[speciei], pos, rho.name())
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
                interpolate(pD[speciei], neg, rho.name())
            )
        );
    }
}


surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

volScalarField rPsi(1.0/psi);
surfaceScalarField rPsi_pos(interpolate(rPsi, pos, "T"));
surfaceScalarField rPsi_neg(interpolate(rPsi, neg, "T"));

surfaceScalarField e_pos(interpolate(e, pos, "T"));
surfaceScalarField e_neg(interpolate(e, neg, "T"));

PtrList<surfaceScalarField> evk_pos(noSolvedVibEq);
PtrList<surfaceScalarField> evk_neg(noSolvedVibEq);

if (downgradeSingleT)
{
    evk_pos.clear();
    evk_neg.clear();
}
else if (downgradeSingleTv)
{
    evk_pos.resize(1);
    evk_neg.resize(1);

    evk_pos.set
    (
        0,
        new surfaceScalarField
        (
            "ev_pos",
            interpolate(ev, pos, "T")
        )
    );

    evk_neg.set
    (
        0,
        new surfaceScalarField
        (
            "ev_neg",
            interpolate(ev, neg, "T")
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
                interpolate(hevk[speciei], pos, "T")
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
                interpolate(hevk[speciei], neg, "T")
            )
        );
    }
}

volScalarField gammatr("gammatr", thermo.Cp_t()/thermo.Cv_t()); // NEW VINCENT (also useful for bdry cdts)
volScalarField c(sqrt(gammatr*rPsi)); // See Candler 2009
Mach = mag(U)/c;

surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

surfaceScalarField cSf_pos //- Intermediate quantity to appear in Eq.8a
(
    "cSf_pos",
    interpolate(c, pos, "T")*mesh.magSf()
);
cSf_pos.setOriented();

surfaceScalarField cSf_neg //- Intermediate quantity to appear in Eq.8b
(
    "cSf_neg",
    interpolate(c, neg, "T")*mesh.magSf()
);
cSf_neg.setOriented();
