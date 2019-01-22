
//------------------------------------------------------------------------------
// Code adapted from Praveen Chandrashekar's solver gFoam
// Personal website: http://cpraveen.github.io/index.html
// gFoam website: https://sourceforge.net/u/cpraveen/
//------------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Lax-Fridrichs flux function
// -----------------------------------------------------------------------------
void LxFFlux 
(
    const bool& downgradeSingleT,
    const bool& downgradeSingleTv,
    const vector& Sf,
    const scalar& cSf_pos,
    const scalar& cSf_neg,
    const scalar& rho_pos,
    const scalar& rho_neg,
    const PtrList<scalar>& pD_pos,
    const PtrList<scalar>& pD_neg,
    const vector& U_pos,
    const vector& U_neg,
    const vector& rhoU_pos,
    const vector& rhoU_neg,
    const scalar& p_pos,
    const scalar& p_neg,
    const PtrList<scalar>& evk_pos,
    const PtrList<scalar>& evk_neg,
    const scalar& e_pos,
    const scalar& e_neg,
    scalar& amaxSf,
    vector& aByU,
    scalar& rhoFlux,
    vector& rhoUFlux,
    PtrList<scalar>& rhoEvkFlux,
    scalar& rhoEFlux
)
{
    // Maximum eigenvalue
    amaxSf = 0.5 * mag((U_pos + U_neg) & Sf) + 
        0.5 * (cSf_pos + cSf_neg);
        
    aByU = (U_pos + U_neg)/2.0;

    const scalar Un_pos = U_pos & Sf;
    const scalar Un_neg = U_neg & Sf;

    const scalar rhoE_pos = rho_pos * (e_pos + 0.5 * magSqr(U_pos));
    const scalar rhoE_neg = rho_neg * (e_neg + 0.5 * magSqr(U_neg));

    rhoFlux = (rho_pos * Un_pos + rho_neg * Un_neg) 
           - amaxSf * (rho_neg - rho_pos);

    rhoUFlux = (p_pos * Sf + rhoU_pos * Un_pos)
            + (p_neg * Sf + rhoU_neg * Un_neg)
            - amaxSf * (rhoU_neg - rhoU_pos);

    rhoEFlux = (rhoE_pos + p_pos) * Un_pos
            + (rhoE_neg + p_neg) * Un_neg
            - amaxSf * (rhoE_neg - rhoE_pos);

    rhoFlux  *= 0.5;
    rhoUFlux *= 0.5;
    rhoEFlux *= 0.5;
}
