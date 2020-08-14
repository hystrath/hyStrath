
//- Eq.7 applied for Psi = rho
phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

//- Eq.7 applied for Phi = rhoU and extra term for pressure
surfaceVectorField phiUp
(
    (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
  + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
);

//- Eq.7 applied for Phi = rhoE and decomposed extra term for pressure
//  (no omega_f for pressure terms, see def of aphiv_pos)
surfaceScalarField phiEp
(
    "phiEp",
    aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
  + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
  + aSf*p_pos - aSf*p_neg
);

// Make flux for pressure-work absolute
if (mesh.moving())
{
    phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
}

PtrList<surfaceScalarField> phiEvk(noSolvedVibEq);

if (downgradeSingleT)
{
    phiEvk.clear();
}
else if (downgradeSingleTv)
{
    phiEvk.resize(1);

    phiEvk.set
    (
        0,
        new surfaceScalarField
        (
            "phiEv",
            aphiv_pos*rho_pos*evk_pos[0]
          + aphiv_neg*rho_neg*evk_neg[0]
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
                "phiEvk_" + Y[speciei].name(),
                aphiv_pos*(pD_pos[speciei]*evk_pos[speciei])
              + aphiv_neg*(pD_neg[speciei]*evk_neg[speciei])
            )
        );
    }
}

