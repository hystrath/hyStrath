surfaceScalarField splitMachDeg1_pos
(
    const surfaceScalarField& interpolatedMach
)
{
    return 0.5*(interpolatedMach+mag(interpolatedMach));
}


surfaceScalarField splitMachDeg1_neg
(
    const surfaceScalarField& interpolatedMach
)
{
    return 0.5*(interpolatedMach-mag(interpolatedMach));
}


surfaceScalarField splitMachDeg2_pos
(
    const surfaceScalarField& interpolatedMach
)
{
    return 0.25*sqr(interpolatedMach+1.0);
}


surfaceScalarField splitMachDeg2_neg
(
    const surfaceScalarField& interpolatedMach
)
{
    return -0.25*sqr(interpolatedMach-1.0);
}


surfaceScalarField splitMachDeg4_pos
(
    const surfaceScalarField& interpolatedMach,
    const scalar AUSM_beta
)
{
    surfaceScalarField sMD4_pos = splitMachDeg1_pos(interpolatedMach);

    forAll(sMD4_pos, intFacei)
    {
        if(mag(interpolatedMach[intFacei]) < 1.0)
        {
            sMD4_pos[intFacei] = splitMachDeg2_pos(interpolatedMach)[intFacei]
                *(1.0-16.0*AUSM_beta*splitMachDeg2_neg(interpolatedMach)[intFacei]);
        }
    }

    forAll(sMD4_pos.boundaryField(), patchi)
    {
        forAll(sMD4_pos.boundaryField()[patchi], facei)
        {
            if(mag(interpolatedMach.boundaryField()[patchi][facei]) < 1.0)
            {
                sMD4_pos.boundaryFieldRef()[patchi][facei] = splitMachDeg2_pos(interpolatedMach).boundaryField()[patchi][facei]
                    *(1.0-16.0*AUSM_beta*splitMachDeg2_neg(interpolatedMach).boundaryField()[patchi][facei]);
            }
        }
    }

    return sMD4_pos;
}


surfaceScalarField splitMachDeg4_neg
(
    const surfaceScalarField& interpolatedMach,
    const scalar AUSM_beta
)
{
    surfaceScalarField sMD4_neg = splitMachDeg1_neg(interpolatedMach);
    
    forAll(sMD4_neg, intFacei)
    {
        if(mag(interpolatedMach[intFacei]) < 1.0)
        {
            sMD4_neg[intFacei] = splitMachDeg2_neg(interpolatedMach)[intFacei]
                *(1.0 + 16.0*AUSM_beta*splitMachDeg2_pos(interpolatedMach)[intFacei]);
        }
    }

    forAll(sMD4_neg.boundaryField(), patchi)
    {
        forAll(sMD4_neg.boundaryField()[patchi], facei)
        {
            if(mag(interpolatedMach.boundaryField()[patchi][facei]) < 1.0)
            {
                sMD4_neg.boundaryFieldRef()[patchi][facei] = splitMachDeg2_neg(interpolatedMach).boundaryField()[patchi][facei]
                    *(1.0 + 16.0*AUSM_beta*splitMachDeg2_pos(interpolatedMach).boundaryField()[patchi][facei]);
            }
        }
    }
    
    return sMD4_neg;
}


surfaceScalarField splitPressureDeg5_pos
(
    const surfaceScalarField& interpolatedMach,
    const surfaceScalarField& AUSM_alpha
)
{
    surfaceScalarField sPD5_pos = splitMachDeg1_pos(interpolatedMach)/(interpolatedMach+1e-06); //- 1e-6: safe-guard for fluid initially at rest

    forAll(sPD5_pos, intFacei)
    {
        if(mag(interpolatedMach[intFacei]) < 1.0)
        {
            sPD5_pos[intFacei] = splitMachDeg2_pos(interpolatedMach)[intFacei]*(2.0-interpolatedMach[intFacei]) 
                - 16.0*AUSM_alpha[intFacei]*interpolatedMach[intFacei]*splitMachDeg2_neg(interpolatedMach)[intFacei];
        }
    }

    forAll(sPD5_pos.boundaryField(), patchi)
    {
        forAll(sPD5_pos.boundaryField()[patchi], facei)
        {
            if(mag(interpolatedMach.boundaryField()[patchi][facei]) < 1.0)
            {
                sPD5_pos.boundaryFieldRef()[patchi][facei] = splitMachDeg2_pos(interpolatedMach).boundaryField()[patchi][facei]*(2.0-interpolatedMach.boundaryField()[patchi][facei]) 
                    - 16.0*AUSM_alpha.boundaryField()[patchi][facei]*interpolatedMach.boundaryField()[patchi][facei]*splitMachDeg2_neg(interpolatedMach).boundaryField()[patchi][facei];
            }
        }
    }
    
    return sPD5_pos;
}


surfaceScalarField splitPressureDeg5_neg
(
    const surfaceScalarField& interpolatedMach,
    const surfaceScalarField& AUSM_alpha
)
{
    surfaceScalarField sPD5_neg = splitMachDeg1_neg(interpolatedMach)/(interpolatedMach+1e-06); //- 1e-6: safe-guard for fluid initially at rest
    
    forAll(sPD5_neg, intFacei)
    {
        if(mag(interpolatedMach[intFacei]) < 1.0)
        {
            sPD5_neg[intFacei] = splitMachDeg2_neg(interpolatedMach)[intFacei]*(-2.0-interpolatedMach[intFacei]) 
                + 16.0*AUSM_alpha[intFacei]*interpolatedMach[intFacei]*splitMachDeg2_pos(interpolatedMach)[intFacei];
        }
    }
    
    forAll(sPD5_neg.boundaryField(), patchi)
    {
        forAll(sPD5_neg.boundaryField()[patchi], facei)
        {
            if(mag(interpolatedMach.boundaryField()[patchi][facei]) < 1.0)
            {
                sPD5_neg.boundaryFieldRef()[patchi][facei] = splitMachDeg2_neg(interpolatedMach).boundaryField()[patchi][facei]*(-2.0-interpolatedMach.boundaryField()[patchi][facei]) 
                    + 16.0*AUSM_alpha.boundaryField()[patchi][facei]*interpolatedMach.boundaryField()[patchi][facei]*splitMachDeg2_pos(interpolatedMach).boundaryField()[patchi][facei];
            }
        }
    }

    return sPD5_neg;
}
