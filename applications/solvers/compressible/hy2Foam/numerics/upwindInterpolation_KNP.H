
surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf()); //- Eq in the text after Eq.6, + direc
surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf()); //- Eq in the text after Eq.6, - direc

surfaceScalarField cSf_pos //- Intermediate quantity to appear in Eq.8a
(
    "cSf_pos",
    fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf()
);
surfaceScalarField cSf_neg //- Intermediate quantity to appear in Eq.8b
(
    "cSf_neg",
    fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf()
);

surfaceScalarField ap //- Eq.8a
(
    "ap",
    max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
);
surfaceScalarField am //- Eq.8b
(
    "am",
    min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
);

surfaceScalarField a_pos("a_pos", ap/(ap - am)); //- Eq.9 KNP

surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

surfaceScalarField aSf("aSf", am*a_pos);

if (fluxScheme == "Tadmor")
{
    aSf = -0.5*amaxSf; //- Eq.10 KT
    a_pos = 0.5; //- Eq.9 KT
}

surfaceScalarField a_neg("a_neg", 1.0 - a_pos); //- (1-alpha) with alpha as in Eq.9 KNP

phiv_pos *= a_pos; //- First factor in the RHS of Eq.7
phiv_neg *= a_neg; //- Second factor in the RHS of Eq.7

surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf); //- Factor of Psi_f+ in factorised Eq.7
surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf); //- Factor of Psi_f- in factorised Eq.7

// Reuse amaxSf for the maximum positive and negative fluxes
// estimated by the central scheme
amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

