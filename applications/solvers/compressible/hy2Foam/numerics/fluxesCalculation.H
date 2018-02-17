
//#include "numerics/upwindInterpolation_KNP.H"
//#include "numerics/upwindInterpolation_AUSM2.H"

// --- Fluxes calculation
//#include "numerics/fluxesCalculation-KNP.H"

if
(
    fluxScheme == "Kurganov" or fluxScheme == "Tadmor" or fluxScheme == "ROE"
 or fluxScheme == "KFVS" or fluxScheme == "LxF"
)
{
    scalar pDListSizes = phiEvk.size();
    if(downgradeSingleT or downgradeSingleTv)
    {
        pDListSizes = 0;
    }

    //- Internal faces
    forAll(mesh.Sf(), facei)
    {
        PtrList<scalar> pDens_pos(pDListSizes);
        PtrList<scalar> pDens_neg(pDListSizes);
        PtrList<scalar> ev_pos(phiEvk.size());
        PtrList<scalar> ev_neg(phiEvk.size());
        PtrList<scalar> phiEv(phiEvk.size());
        
        forAll(pDens_pos, speciei)
        { 
            pDens_pos.set
            (
                speciei,
                new scalar(pD_pos[speciei][facei])
            );
            
            pDens_neg.set
            (
                speciei,
                new scalar(pD_neg[speciei][facei])
            );
        }
        
        forAll(ev_pos, speciei)
        {
            ev_pos.set
            (
                speciei,
                new scalar(evk_pos[speciei][facei])
            );
            
            ev_neg.set
            (
                speciei,
                new scalar(evk_neg[speciei][facei])
            );
            
            phiEv.set
            (
                speciei,
                new scalar(phiEvk[speciei][facei])
            );
        }
        
        fluxMap[fluxScheme]
        (
            downgradeSingleT,
            downgradeSingleTv,
            mesh.Sf()[facei],
            cSf_pos[facei],
            cSf_neg[facei],
            rho_pos[facei],
            rho_neg[facei],
            pDens_pos,
            pDens_neg,
            U_pos[facei],
            U_neg[facei],
            rhoU_pos[facei],
            rhoU_neg[facei],
            p_pos[facei],
            p_neg[facei],
            ev_pos,
            ev_neg,
            e_pos[facei],
            e_neg[facei],
            amaxSf[facei],
            aByU[facei],
            phi[facei],
            phiUp[facei],
            phiEv,
            phiEp[facei]
        );
        
        forAll(phiEvk, speciei)
        {
            phiEvk[speciei][facei] = phiEv[speciei];
        }
    }
    
    //- Boundary faces
    forAll(mesh.boundaryMesh(), patchi)
    {
        if(mesh.boundaryMesh()[patchi].type() != "empty")
        {
            forAll(mesh.boundaryMesh()[patchi], facei)
            {
                PtrList<scalar> ppDens_pos(pDListSizes);
                PtrList<scalar> ppDens_neg(pDListSizes);
                PtrList<scalar> pev_pos(phiEvk.size());
                PtrList<scalar> pev_neg(phiEvk.size());
                PtrList<scalar> pphiEv(phiEvk.size());
                
                forAll(ppDens_pos, speciei)
                { 
                    ppDens_pos.set
                    (
                        speciei,
                        new scalar
                        (
                            pD_pos[speciei].boundaryField()[patchi][facei]
                        )
                    );
                    
                    ppDens_neg.set
                    (
                        speciei,
                        new scalar
                        (
                            pD_neg[speciei].boundaryField()[patchi][facei]
                        )
                    );
                }
                
                forAll(pev_pos, speciei)
                {
                    pev_pos.set
                    (
                        speciei,
                        new scalar
                        (
                            evk_pos[speciei].boundaryField()[patchi][facei]
                        )
                    );
                    
                    pev_neg.set
                    (
                        speciei,
                        new scalar
                        (
                            evk_neg[speciei].boundaryField()[patchi][facei]
                        )
                    );
                    
                    pphiEv.set
                    (
                        speciei,
                        new scalar
                        (
                            phiEvk[speciei].boundaryField()[patchi][facei]
                        )
                    );
                }
        
                fluxMap[fluxScheme]
                (
                    downgradeSingleT,
                    downgradeSingleTv,
                    mesh.boundaryMesh()[patchi].faceAreas()[facei],
                    cSf_pos.boundaryField()[patchi][facei],
                    cSf_neg.boundaryField()[patchi][facei],
                    rho_pos.boundaryField()[patchi][facei],
                    rho_neg.boundaryField()[patchi][facei],
                    ppDens_pos,
                    ppDens_neg,
                    U_pos.boundaryField()[patchi][facei],
                    U_neg.boundaryField()[patchi][facei],
                    rhoU_pos.boundaryField()[patchi][facei],
                    rhoU_neg.boundaryField()[patchi][facei],
                    p_pos.boundaryField()[patchi][facei],
                    p_neg.boundaryField()[patchi][facei],
                    pev_pos,
                    pev_neg,
                    e_pos.boundaryField()[patchi][facei],
                    e_neg.boundaryField()[patchi][facei],
                    amaxSf.boundaryFieldRef()[patchi][facei],
                    aByU.boundaryFieldRef()[patchi][facei],
                    phi.boundaryFieldRef()[patchi][facei],
                    phiUp.boundaryFieldRef()[patchi][facei],
                    pphiEv,
                    phiEp.boundaryFieldRef()[patchi][facei]
                );
                
                forAll(phiEvk, speciei)
                {
                    phiEvk[speciei].boundaryFieldRef()[patchi][facei] =
                        pphiEv[speciei];
                }
            }
        }
    }
}

