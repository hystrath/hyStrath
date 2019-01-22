
/*surfaceScalarField phiv_pos("phiv_pos", (U_pos & mesh.Sf())/mesh.magSf()); //- Step II, + direction
surfaceScalarField phiv_neg("phiv_neg", (U_neg & mesh.Sf())/mesh.magSf()); //- Step II, - direction*/

surfaceScalarField phiv_pos("phiv_pos", (U_pos & mesh.Sf())); //- Step II, + direction
surfaceScalarField phiv_neg("phiv_neg", (U_neg & mesh.Sf())); //- Step II, - direction

/*Info << "phiv_pos2: " << phiv_pos[2] << tab << "phiv_neg" << phiv_neg[2] << endl;
Info << "phiv_pos3: " << phiv_pos[3] << tab << "phiv_neg" << phiv_neg[3] << endl;
Info << "phiv_pos4: " << phiv_pos[4] << tab << "phiv_neg" << phiv_neg[4] << endl;
Info << "phiv_pos5: " << phiv_pos[5] << tab << "phiv_neg" << phiv_neg[5] << endl;*/

surfaceScalarField c_pos //- Step I, + direction
(
    "c_pos",
    fvc::interpolate(c, pos, "reconstruct(T)")
);

surfaceScalarField c_neg //- Step I, - direction
(
    "c_neg",
    fvc::interpolate(c, neg, "reconstruct(T)")
);

surfaceScalarField cf //- Step I, cf definition no 1
(
    "cf",
    0.5*(c_pos + c_neg)
);

Info << "cf dim: " << cf.dimensions() << endl;
Info << "phi dim: " << phi.dimensions() << endl;

surfaceScalarField Mach_pos //- Eq.69
(
    "Mach_pos",
    phiv_pos/cf
);

surfaceScalarField Mach_neg
(
    "Mach_neg",
    phiv_neg/cf
);

Info << "Mach_pos: " << Mach_pos.dimensions() << tab << Mach_pos[2] << tab << "Mach_neg: " << Mach_neg[2] << endl;
Info << "Mach_pos: " << Mach_pos.boundaryField() << endl;

surfaceScalarField meanSqrMach //- Step III
(
    "meanSqrMach",
    (sqr(phiv_pos) + sqr(phiv_neg))/(2.0*sqr(cf))
);

Info << "meanSqrMach: " << meanSqrMach[1] << tab << meanSqrMach[2] << endl;

/*surfaceScalarField meanSqrMachZero //- Step IV
(
    "meanSqrMachZero",
    min(1.0, max(max(meanSqrMach, sqrMachInf), 1e-6)) //- 1e-6: safe-guard for fluid initially at rest
);

Info << "meanSqrMachZero: " << meanSqrMachZero[1] << tab << meanSqrMachZero[2] << endl;*/

surfaceScalarField scalingFunction //- Step V
(
    "fa",
    pos
    //meanSqrMachZero*(2.0-meanSqrMachZero)
);

/*surfaceScalarField scalingFunction //- Step V
(
    "fa",
    meanSqrMachZero*(2.0-meanSqrMachZero)
);


Info << "scalingFunction: " << scalingFunction[1] << tab << scalingFunction[2] << endl;*/

scalar AUSM_Kp = 0.25;
scalar AUSM_Ku = 0.75;
scalar AUSM_sigma = 1.0;
scalar AUSM_beta = 1.0/8.0;

surfaceScalarField AUSM_alpha
(
    "AUSM_alpha",
    3.0/16.0*(-4.0+5.0*sqr(scalingFunction))
);

Info << "Mach_pos boundary: " << Mach_pos.boundaryField() << endl;
Info << "Mach_neg boundary: " << Mach_neg.boundaryField() << endl;

surfaceScalarField Machf //- Step VI
(
    "Machf",
    //splitMachDeg4_pos(Mach_neg, AUSM_beta) + splitMachDeg4_neg(Mach_pos, AUSM_beta)
    splitMachDeg4_neg(Mach_neg, AUSM_beta) + splitMachDeg4_pos(Mach_pos, AUSM_beta)
);

/*Machf -= AUSM_Kp/scalingFunction*max(1.0-AUSM_sigma*meanSqrMach, 0.0)*(p_pos-p_neg)
  / (0.5*(rho_pos+rho_neg)*sqr(cf));*/
  
Info << "Machf" << tab << Machf[2] << tab << Machf[3] << endl;

surfaceScalarField pf 
(
    "pf",
    //splitPressureDeg5_pos(Mach_neg, AUSM_alpha)*p_neg + splitPressureDeg5_neg(Mach_pos, AUSM_alpha)*p_pos
    splitPressureDeg5_neg(Mach_neg, AUSM_alpha)*p_neg + splitPressureDeg5_pos(Mach_pos, AUSM_alpha)*p_pos
);

/*pf -= AUSM_Ku*(rho_pos+rho_neg)*(phiv_pos-phiv_neg)*(scalingFunction*cf)
   *splitPressureDeg5_pos(Mach_neg, AUSM_alpha)*splitPressureDeg5_neg(Mach_pos, AUSM_alpha);*/

Info << "hello9" << endl;

// Use amaxSf for the maximum positive and negative fluxes
// estimated by the AUSM+-up scheme
surfaceScalarField amaxSf = max(mag(phiv_pos), mag(phiv_neg))*mesh.magSf();

Info << "hello10" << endl;

