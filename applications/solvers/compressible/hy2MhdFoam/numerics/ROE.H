#include <numeric>

//------------------------------------------------------------------------------
// Code adapted from Praveen Chandrashekar's solver gFoam
// Personal website: http://cpraveen.github.io/index.html
// gFoam website: https://sourceforge.net/u/cpraveen/
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// ROE flux function
// There is no entropy fix
//------------------------------------------------------------------------------
void ROEFlux 
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
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf)
        + 0.5*(cSf_pos + cSf_neg);
        
    aByU = (U_pos + U_neg)/2.0;    
        
    const scalar area = mag(Sf);
    const vector unit_normal = Sf/area;
    
    const scalar sum_evk_pos = 
        std::accumulate(evk_pos.begin(), evk_pos.end(), 0.0);
        
    const scalar etr_pos = e_pos - sum_evk_pos;

    const scalar gammaMinusOne = p_pos/(rho_pos*etr_pos);

    //- Enthalpy
    const scalar h_pos  = e_pos + p_pos/rho_pos + 0.5 * magSqr(U_pos);
    const scalar h_neg  = e_neg + p_neg/rho_neg + 0.5 * magSqr(U_neg);

    const scalar rho_pos_sqrt = Foam::sqrt(rho_pos);
    const scalar rho_neg_sqrt = Foam::sqrt(rho_neg);
    const scalar fact_pos = rho_pos_sqrt/(rho_pos_sqrt + rho_neg_sqrt);
    const scalar fact_neg = 1.0 - fact_pos;

    //- Roe average state
    const scalar density  = rho_pos_sqrt*rho_neg_sqrt;
    const vector velocity = U_pos*fact_pos + U_neg*fact_neg;
    const scalar h = h_pos*fact_pos + h_neg*fact_neg;
 
    const scalar vel_normal = velocity & unit_normal;
    const scalar c = Foam::sqrt(gammaMinusOne*(h - 0.5*magSqr(velocity))); // TODO

    const scalar dp = p_neg - p_pos;
    const scalar vel_pos_normal = U_pos & unit_normal;
    const scalar vel_neg_normal = U_neg & unit_normal;
    const scalar dV = vel_neg_normal - vel_pos_normal;
    
    if(vel_normal > 0.0)
    {
        const scalar lambda = vel_normal - c;
        const scalar coeff  = 0.5*(dp - density*c*dV)/Foam::sqr(c);
        const scalar factor = Foam::min(lambda, 0.0)*coeff;

        //- Left flux
        rhoFlux  = rho_pos*vel_pos_normal;
        rhoUFlux = unit_normal*p_pos + rho_pos*U_pos*vel_pos_normal;
        rhoEFlux = h_pos*rhoFlux;

        //- Upwind term
        rhoFlux  += factor;
        rhoUFlux += (velocity - unit_normal*c)*factor;
        rhoEFlux += (h - c*vel_normal)*factor;
    }
    else
    {
        const scalar lambda = vel_normal + c;
        const scalar coeff  = 0.5*(dp + density*c*dV)/Foam::sqr(c);
        const scalar factor = Foam::max(lambda, 0.0)*coeff;

        //- Right flux
        rhoFlux  = rho_neg*vel_neg_normal;
        rhoUFlux = unit_normal*p_neg + rho_neg*U_neg*vel_neg_normal;
        rhoEFlux = h_neg*rhoFlux;

        //- Upwind term
        rhoFlux  -= factor;
        rhoUFlux -= (velocity + unit_normal*c)*factor;
        rhoEFlux -= (h + c*vel_normal)*factor;
    }

    rhoFlux  *= area;
    rhoUFlux *= area;
    rhoEFlux *= area;
}
