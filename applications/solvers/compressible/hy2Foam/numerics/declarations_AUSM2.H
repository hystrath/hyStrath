
phi = cf*(max(Machf, vvvv_zero)*rho_pos + min(Machf, vvvv_zero)*rho_neg); //- Eq.74

Info << "phi: " << phi << endl;

surfaceVectorField phiUp
(
    max(phi, vvv_zero)*U_pos + min(phi, vvv_zero)*U_neg
);

Info << "phiUp without p: " << phiUp << endl;

phiUp += pf*mesh.Sf();

Info << "phiUp with p: " << phiUp << endl;

surfaceScalarField phiEp
(
    "phiEp",
    max(phi, vvv_zero)*((e_pos + 0.5*magSqr(U_pos)) + p_pos/rho_pos)
  + min(phi, vvv_zero)*((e_neg + 0.5*magSqr(U_neg)) + p_neg/rho_neg)

    /*
    aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
  + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
  + aSf*p_pos - aSf*p_neg*/
);

// Make flux for pressure-work absolute
if (mesh.moving())
{
    phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
}

Info << "phiEp: " << phiEp << endl;

// NEW VINCENT ****************************************************************
PtrList<surfaceScalarField> phiEvk(noSolvedVibEq);

if(downgradeSingleT) // NEW VINCENT 05-16/08/2016
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
            "phiEv",
            max(phi, vvv_zero)*rho_pos*evk_pos[0]
          + min(phi, vvv_zero)*rho_neg*evk_neg[0]
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
                max(phi, vvv_zero)*(pD_pos[speciei]*evk_pos[speciei])
              + min(phi, vvv_zero)*(pD_neg[speciei]*evk_neg[speciei])
            )
        );
    }
}

// END NEW VINCENT ************************************************************

volScalarField muEff(turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

if(hyDev)
{
    Zrot = ZrotInf/(1.0 + sqr(pi)/2.0*sqrt(TZ/thermo.T()) + (sqr(pi)/4.0 + pi)*TZ/thermo.T());
    Z = 3.0*pi/(4.0*(thermo.zetar()+3.0))*Zrot;

    mub = mu*(2.0*thermo.zetar()*Z)/(3.0*(thermo.zetar()+3.0));
    bulkToShearRatio = mub/mu;

    tauMC = tauMC + mub*tr(fvc::grad(U))*I;
}

