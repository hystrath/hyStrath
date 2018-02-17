
//------------------------------------------------------------------------------
// Code adapted from Praveen Chandrashekar's solver gFoam
// Personal website: http://cpraveen.github.io/index.html
// gFoam website: https://sourceforge.net/u/cpraveen/
//------------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Kinetic split fluxes
// sign = +1 gives positive flux
// sign = -1 gives negative flux
// -----------------------------------------------------------------------------
void KFVSSplitFlux 
(
    const scalar& sign,
    const vector& normal,
    const scalar& rho,
    const vector& U,
    const scalar& p,
    const scalar& e,
    scalar& rhoFlux,  
    vector& rhoUFlux, 
    scalar& rhoEFlux
)
{
   const scalar vn   = (U & normal);
   const scalar beta = 0.5*rho/p;
   const scalar s    = vn*Foam::sqrt(beta);
   const scalar A    = 0.5*(1.0 + sign*Foam::erf(s));
   const scalar B    = 0.5*sign*Foam::exp(-s*s) / 
                       Foam::sqrt(constant::mathematical::pi*beta);
   const scalar ufact= vn*A + B;
   const scalar rhoE = rho*(e + 0.5*magSqr(U));

   rhoFlux  = rho*ufact;
   rhoUFlux = p*A*normal + rho*U*ufact;
   rhoEFlux = (rhoE + p)*vn*A + (rhoE + 0.5*p)*B;
}

// -----------------------------------------------------------------------------
// Kinetic flux function
// -----------------------------------------------------------------------------
void KFVSFlux
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
    
    scalar rhoFlux_pos, rhoFlux_neg;
    vector rhoUFlux_pos, rhoUFlux_neg;
    scalar rhoEFlux_pos, rhoEFlux_neg;

    //- unit normal
    const scalar area = mag(Sf);
    const vector normal = Sf/area;

    //- positive flux
    KFVSSplitFlux
    (
        1.0,
        normal,
        rho_pos,
        U_pos,
        p_pos,
        e_pos,
        rhoFlux_pos,  
        rhoUFlux_pos, 
        rhoEFlux_pos
    );

    //- negative flux
    KFVSSplitFlux
    (
       -1.0,
        normal,
        rho_neg,
        U_neg,
        p_neg,
        e_neg,
        rhoFlux_neg,  
        rhoUFlux_neg, 
        rhoEFlux_neg
    );

    rhoFlux  = (rhoFlux_pos + rhoFlux_neg)*area;
    rhoUFlux = (rhoUFlux_pos + rhoUFlux_neg)*area;
    rhoEFlux = (rhoEFlux_pos + rhoEFlux_neg)*area;
}
