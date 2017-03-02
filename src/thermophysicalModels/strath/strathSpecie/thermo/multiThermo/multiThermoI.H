/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type>::multiThermo
(
    const MultiThermo& sp
)
:
    MultiThermo(sp)
{}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (multiThermo<MultiThermo, Type>::*F)(const scalar, const scalar) const,
    scalar (multiThermo<MultiThermo, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (multiThermo<MultiThermo, Type>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "multiThermo<MultiThermo, Type>::T(scalar f, scalar T0, "
                "scalar (multiThermo<MultiThermo, Type>::*F)"
                "(const scalar) const, "
                "scalar (multiThermo<MultiThermo, Type>::*dFdT)"
                "(const scalar) const, "
                "scalar (multiThermo<MultiThermo, Type>::*limit)"
                "(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::T
(
    label vibMode,
    scalar f,
    scalar p,
    scalar T0,
    scalar (multiThermo<MultiThermo, Type>::*F)(const label, const scalar, const scalar) const,
    scalar (multiThermo<MultiThermo, Type>::*dFdT)(const label, const scalar, const scalar)
        const,
    scalar (multiThermo<MultiThermo, Type>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(vibMode, p, Test) - f)/(this->*dFdT)(vibMode, p, Test));   

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "multiThermo<MultiThermo, Type>::T(label vibMode, scalar f, scalar T0, "
                "scalar (multiThermo<MultiThermo, Type>::*F)"
                "(const label, const scalar) const, "
                "scalar (multiThermo<MultiThermo, Type>::*dFdT)"
                "(const label, const scalar) const, "
                "scalar (multiThermo<MultiThermo, Type>::*limit)"
                "(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type>::multiThermo
(
    const word& name,
    const multiThermo& st
)
:
    MultiThermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
inline Foam::word
Foam::species::multiThermo<MultiThermo, Type>::heName()
{
    return Type<multiThermo<MultiThermo, Type> >::name();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::he(const scalar p, const scalar Tt) const
{
    return Type<multiThermo<MultiThermo, Type> >::he(*this, p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::he(const scalar p, const scalar Tt, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::he(*this, p, Tt, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::het(const scalar p, const scalar Tt) const
{
    return Type<multiThermo<MultiThermo, Type> >::het(*this, p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hev_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::hev_mode(*this, vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hev(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::hev(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::heel(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::heel(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hevel_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::hevel_mode(*this, vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hevel(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::hevel(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cv(const scalar p, const scalar Tt) const
{
    return this->cv_t(p, Tt) + this->cv_v(p, Tt) + this->cv_el(p, Tt) + this->cv_e(p, Tt); // NEW CORRECTION VINCENT 21/01/2016
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cv(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->cv_t(p, Tt) + this->cv_v(p, Tv) + this->cv_el(p, Tv) + this->cv_e(p, Tv); // NEW CORRECTION VINCENT 21/01/2016
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cv_vel(const scalar p, const scalar Tv) const
{
    return this->cv_v(p, Tv) + this->cv_el(p, Tv) + this->cv_e(p, Tv); // NEW CORRECTION VINCENT 21/01/2016
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cv_vel_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    if(vibMode == 0) // NEW VINCENT 12/03/2016
    {
        return this->cv_v_mode(vibMode, p, Tv) + this->cv_el(p, Tv) + this->cv_e(p, Tv);
    }
    else
    {
        return this->cv_v_mode(vibMode, p, Tv);
    }
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp(const scalar p, const scalar Tt) const
{
    return this->cv(p, Tt) + this->cpMcv(p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->cv(p, Tt, Tv) + this->cpMcv(p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp_t(const scalar p, const scalar Tt) const
{
    return this->cv_t(p, Tt) + this->cpMcv(p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp_v_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->cv_v_mode(vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp_v(const scalar p, const scalar Tv) const
{
    return this->cv_v(p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp_el(const scalar p, const scalar Tv) const
{
    return this->cv_el(p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cp_vel(const scalar p, const scalar Tv) const
{
    return this->cv_vel(p, Tv);
}
// END NEW VINCENT ************************************************************


/*template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::cpv(const scalar p, const scalar Tt, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::cpv(*this, p, Tt, Tv);
}*/


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::gamma(const scalar p, const scalar Tt, const scalar Tv) const
{
    scalar cv = this->cv(p, Tt, Tv);
    return (cv + this->cpMcv(p, Tt))/cv;
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::es(const scalar p, const scalar Tt) const
{
    return this->ets(p, Tt) + this->evs(p, Tt) + this->eels(p, Tt) + this->ees(p, Tt); // NEW CORRECTION VINCENT 21/01/2016
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::es(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->ets(p, Tt) + this->evs(p, Tv) + this->eels(p, Tv) + this->ees(p, Tv); // NEW CORRECTION VINCENT 21/01/2016
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::evels_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    if(vibMode == 0)
    {
        return this->evs_mode(vibMode, p, Tv) + this->eels(p, Tv) + this->ees(p, Tv);
    }
    else
    {
        return this->evs_mode(vibMode, p, Tv);
    }
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::evels(const scalar p, const scalar Tv) const
{
    return this->evs(p, Tv) + this->eels(p, Tv) + this->ees(p, Tv); // NEW CORRECTION VINCENT 21/01/2016
}


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hs(const scalar p, const scalar Tt) const
{
    return this->hts(p, Tt) + this->hvs(p, Tt) + this->hels(p, Tt) + this->ees(p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hs(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->hts(p, Tt) + this->hvs(p, Tv) + this->hels(p, Tv) + this->ees(p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hts(const scalar p, const scalar Tt) const
{
    return this->ets(p, Tt) + p*this->W()/this->rho(p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hvs_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->evs_mode(vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hvs(const scalar p, const scalar Tv) const
{
    return this->evs(p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hels(const scalar p, const scalar Tv) const
{
    return this->eels(p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hvels_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->evels_mode(vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::hvels(const scalar p, const scalar Tv) const
{
    return this->evels(p, Tv);
}
// END NEW VINCENT ************************************************************


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::ea(const scalar p, const scalar Tt, const scalar Tv) const
{
    return  this->es(p, Tt, Tv) + this->hc();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::g(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->ha(p, Tt, Tv) - Tt*this->s(p, Tt, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::a(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->ea(p, Tt, Tv) - Tt*this->s(p, Tt, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cpv(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->cpv(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->cp(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv(const scalar p, const scalar Tt) const
{
    return this->cv(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->cv(p, Tt, Tv)/this->W();
}


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_trans(const scalar p, const scalar Tt) const
{
    return this->cv_trans(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_r(const scalar p, const scalar Tt) const
{
    return this->cv_r(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_t(const scalar p, const scalar Tt) const
{
    return this->cv_t(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_v_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->cv_v_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_v(const scalar p, const scalar Tv) const
{
    return this->cv_v(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_el(const scalar p, const scalar Tv) const
{
    return this->cv_el(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_vel_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->cv_vel_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cv_vel(const scalar p, const scalar Tv) const
{
    return this->cv_vel(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_t(const scalar p, const scalar Tt) const
{
    return this->cp_t(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_v_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->cp_v_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_v(const scalar p, const scalar Tv) const
{
    return this->cp_v(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_el(const scalar p, const scalar Tv) const
{
    return this->cp_el(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_vel(const scalar p, const scalar Tv) const
{
    return this->cp_vel(p, Tv)/this->W();
}

template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Cp_vel_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->cp_vel_mode(vibMode, p, Tv)/this->W();
}
// END NEW VINCENT ************************************************************


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HE(const scalar p, const scalar Tt) const
{
    return Type<multiThermo<MultiThermo, Type> >::HE(*this, p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HE(const scalar p, const scalar Tt, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HE(*this, p, Tt, Tv);
}


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEt(const scalar p, const scalar Tt) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEt(*this, p, Tt);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEv(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEv(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEv_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEv_mode(*this, vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEel(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEel(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEvel_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEvel_mode(*this, vibMode, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::HEvel(const scalar p, const scalar Tv) const
{
    return Type<multiThermo<MultiThermo, Type> >::HEvel(*this, p, Tv);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hts(const scalar p, const scalar Tt) const
{
    return this->hts(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hvs_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->hvs_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hvs(const scalar p, const scalar Tv) const
{
    return this->hvs(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hels(const scalar p, const scalar Tv) const
{
    return this->hels(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hvels(const scalar p, const scalar Tv) const
{
    return this->hvels(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hvels_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->hvels_mode(vibMode, p, Tv)/this->W();
}
// END NEW VINCENT ************************************************************


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hs(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->hs(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Hc() const
{
    return this->hc()/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Ha(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->ha(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::S(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->s(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::E(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->e(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Es(const scalar p, const scalar Tt) const
{
    return this->es(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Es(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->es(p, Tt, Tv)/this->W();
}

// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Ets(const scalar p, const scalar Tt) const
{
    return this->ets(p, Tt)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Evs_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->evs_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Evs(const scalar p, const scalar Tv) const
{
    return this->evs(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Eels(const scalar p, const scalar Tv) const
{
    return this->eels(p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Evels_mode(const label vibMode, const scalar p, const scalar Tv) const
{
    return this->evels_mode(vibMode, p, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Evels(const scalar p, const scalar Tv) const
{
    return this->evels(p, Tv)/this->W();
}
// END NEW VINCENT ************************************************************


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Ea(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->ea(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::G(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->g(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::A(const scalar p, const scalar Tt, const scalar Tv) const
{
    return this->a(p, Tt, Tv)/this->W();
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->nMoles()*this->g(p, T)/(this->RR*T);

    if (arg < 600.0)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar
Foam::species::multiThermo<MultiThermo, Type>::Kc(const scalar p, const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(this->Pstd/(this->RR*T), this->nMoles());
    }
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(this->Pstd/p, this->nMoles());
    }
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*this->Pstd/p, this->nMoles());
    }
}


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TtHEt
(
    const scalar het,
    const scalar p,
    const scalar T0
) const
{
    return Type<multiThermo<MultiThermo, Type> >::TtHEt(*this, het, p, T0);
}

template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvHEv
(
    const scalar hev,
    const scalar p,
    const scalar T0
) const
{
    return Type<multiThermo<MultiThermo, Type> >::TvHEv(*this, hev, p, T0);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvHEv_mode
(
    const label vibMode,
    const scalar hev,
    const scalar p,
    const scalar T0
) const
{
    return Type<multiThermo<MultiThermo, Type> >::TvHEv_mode(*this, vibMode, hev, p, T0);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvelHEvel
(
    const scalar hevel,
    const scalar p,
    const scalar T0
) const
{
    return Type<multiThermo<MultiThermo, Type> >::TvelHEvel(*this, hevel, p, T0);
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvelHEvel_mode
(
    const label vibMode,
    const scalar hevel,
    const scalar p,
    const scalar T0
) const
{
    return Type<multiThermo<MultiThermo, Type> >::TvelHEvel_mode(*this, vibMode, hevel, p, T0);
}
// END NEW VINCENT ************************************************************


// DELETEION VINCENT **********************************************************
/*template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Hs,
        &multiThermo<MultiThermo, Type>::Cp,
        &multiThermo<MultiThermo, Type>::limit
    );
}*/
// END DELETION VINCENT *******************************************************


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TtHts
(
    const scalar hts,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hts,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Hts,
        &multiThermo<MultiThermo, Type>::Cp_t,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvHvs
(
    const scalar hvs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hvs,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Hvs,
        &multiThermo<MultiThermo, Type>::Cp_v,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvelHvels
(
    const scalar hvels,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hvels,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Hvels,
        &multiThermo<MultiThermo, Type>::Cp_vel,
        &multiThermo<MultiThermo, Type>::limit
    );
}
// END NEW VINCENT ************************************************************


// DELETION VINCENT ***********************************************************
/*template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::THa
(
    const scalar ha,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Ets,
        &multiThermo<MultiThermo, Type>::Cv_t,
        &multiThermo<MultiThermo, Type>::limit
    );
}*/
// END DELETION VINCENT *******************************************************


// DELETION VINCENT ***********************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Es,
        &multiThermo<MultiThermo, Type>::Cv,
        &multiThermo<MultiThermo, Type>::limit
    );
}
// END DELETION VINCENT *******************************************************


// NEW VINCENT ****************************************************************
template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TtEts
(
    const scalar ets,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ets,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Ets,
        &multiThermo<MultiThermo, Type>::Cv_t,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvEvs_mode
(
    const label vibMode,
    const scalar evs_mode,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        vibMode,
        evs_mode,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Evs_mode,
        &multiThermo<MultiThermo, Type>::Cv_v_mode,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvEvs
(
    const scalar evs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        evs,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Evs,
        &multiThermo<MultiThermo, Type>::Cv_v,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvelEvels_mode
(
    const label vibMode,
    const scalar evels_mode,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        vibMode,
        evels_mode,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Evels_mode,
        &multiThermo<MultiThermo, Type>::Cv_vel_mode,
        &multiThermo<MultiThermo, Type>::limit
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TvelEvels
(
    const scalar evels,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        evels,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Evels,
        &multiThermo<MultiThermo, Type>::Cv_vel,
        &multiThermo<MultiThermo, Type>::limit
    );
}
// END NEW VINCENT ************************************************************


// DELETION VINCENT ***********************************************************
/*template<class MultiThermo, template<class> class Type>
inline Foam::scalar Foam::species::multiThermo<MultiThermo, Type>::TEa
(
    const scalar ea,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
        p,
        T0,
        &multiThermo<MultiThermo, Type>::Ea,
        &multiThermo<MultiThermo, Type>::Cv,
        &multiThermo<MultiThermo, Type>::limit
    );
}*/
// END DELETION VINCENT *******************************************************


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
inline void Foam::species::multiThermo<MultiThermo, Type>::operator+=
(
    const multiThermo<MultiThermo, Type>& st
)
{
    MultiThermo::operator+=(st);
}


template<class MultiThermo, template<class> class Type>
inline void Foam::species::multiThermo<MultiThermo, Type>::operator-=
(
    const multiThermo<MultiThermo, Type>& st
)
{
    MultiThermo::operator-=(st);
}


template<class MultiThermo, template<class> class Type>
inline void Foam::species::multiThermo<MultiThermo, Type>::operator*=(const scalar s)
{
    MultiThermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type> Foam::species::operator+
(
    const multiThermo<MultiThermo, Type>& st1,
    const multiThermo<MultiThermo, Type>& st2
)
{
    return multiThermo<MultiThermo, Type>
    (
        static_cast<const MultiThermo&>(st1) + static_cast<const MultiThermo&>(st2)
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type> Foam::species::operator-
(
    const multiThermo<MultiThermo, Type>& st1,
    const multiThermo<MultiThermo, Type>& st2
)
{
    return multiThermo<MultiThermo, Type>
    (
        static_cast<const MultiThermo&>(st1) - static_cast<const MultiThermo&>(st2)
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type> Foam::species::operator*
(
    const scalar s,
    const multiThermo<MultiThermo, Type>& st
)
{
    return multiThermo<MultiThermo, Type>
    (
        s*static_cast<const MultiThermo&>(st)
    );
}


template<class MultiThermo, template<class> class Type>
inline Foam::species::multiThermo<MultiThermo, Type> Foam::species::operator==
(
    const multiThermo<MultiThermo, Type>& st1,
    const multiThermo<MultiThermo, Type>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
