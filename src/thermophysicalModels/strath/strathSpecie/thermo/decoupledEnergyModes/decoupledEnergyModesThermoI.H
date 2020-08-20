/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "decoupledEnergyModesThermo.H"
#include "advancedSpecie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline
Foam::decoupledEnergyModesThermo<EquationOfState>::decoupledEnergyModesThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const typename decoupledEnergyModesThermo<EquationOfState>::cvCoeffsArray& decoupledCvCoeffs,
    const typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray& vibrationalList, // NEW VINCENT 24/02/2016
    const typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray& electronicList
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        decoupledCvCoeffs_[coefLabel] = decoupledCvCoeffs[coefLabel];
    }

    // NEW VINCENT 23/02/2016 *************************************************
    label nCoeffsMin = min(vibrationalList_.size(), vibrationalList.size());
    label nCoeffsMax = max(vibrationalList_.size(), vibrationalList.size());
    if (nCoeffsMin != nCoeffsMax)
    {
        if (nCoeffsMin != 0)
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMin; coefLabel++)
            {
                vibrationalList_[coefLabel] = vibrationalList[coefLabel];
            }
            for (register label coefLabel=nCoeffsMin; coefLabel<nCoeffsMax; coefLabel++)
            {
                vibrationalList_.append(0.0);
            }
        }
        else
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMax; coefLabel++)
            {
                vibrationalList_.append(vibrationalList[coefLabel]);
            }
        }
    }
    // END NEW VINCENT 23/02/2016 *********************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(electronicList_.size(), electronicList.size());
    nCoeffsMax = max(electronicList_.size(), electronicList.size());
    //Info << "1/ "<< "nCoeffsMin: "<< nCoeffsMin<< ", " << "nCoeffsMax: "<< nCoeffsMax << endl;
    //Info << "1/ "<< "electronicList.size(): "<< electronicList.size()<< ", " << "electronicList_.size(): "<< electronicList_.size() << endl;
    if (nCoeffsMin != nCoeffsMax)
    {
        if (nCoeffsMin != 0)
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMin; coefLabel++)
            {
                electronicList_[coefLabel] = electronicList[coefLabel];
            }
            for (register label coefLabel=nCoeffsMin; coefLabel<nCoeffsMax; coefLabel++)
            {
                electronicList_.append(0.0);
            }
        }
        else
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMax; coefLabel++)
            {
                electronicList_.append(electronicList[coefLabel]);
            }
        }
    }
    // END NEW VINCENT 11/02/2016 *********************************************

}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::cvCoeffsArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::cvCoeffs
(
    const scalar T
) const
{
    return decoupledCvCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::dynamicArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::vibrationalCoeffs // NEW VINCENT 24/02/2016
(
    const scalar T
) const
{
    return vibrationalList_;
}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::dynamicArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::electronicCoeffs
(
    const scalar T
) const
{
    return electronicList_;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::decoupledEnergyModesThermo<EquationOfState>::decoupledEnergyModesThermo
(
    const word& name,
    const decoupledEnergyModesThermo& dem
)
:
    EquationOfState(name, dem),
    Tlow_(dem.Tlow_),
    Thigh_(dem.Thigh_)
{
    for (register label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        decoupledCvCoeffs_[coefLabel] = dem.decoupledCvCoeffs_[coefLabel];
    }

    // NEW VINCENT 23/02/2016 *************************************************
    label nCoeffsMin = min(vibrationalList_.size(), dem.vibrationalList_.size());
    label nCoeffsMax = max(vibrationalList_.size(), dem.vibrationalList_.size());
    if (nCoeffsMin != nCoeffsMax)
    {
        if (nCoeffsMin != 0)
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMin; coefLabel++)
            {
                vibrationalList_[coefLabel] = dem.vibrationalList_[coefLabel];
            }
            for (register label coefLabel=nCoeffsMin; coefLabel<nCoeffsMax; coefLabel++)
            {
                vibrationalList_.append(0.0);
            }
        }
        else
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMax; coefLabel++)
            {
                vibrationalList_.append(dem.vibrationalList_[coefLabel]);
            }
        }
    }
    // END NEW VINCENT 23/02/2016 *********************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(electronicList_.size(), dem.electronicList_.size());
    nCoeffsMax = max(electronicList_.size(), dem.electronicList_.size());
    if (nCoeffsMin != nCoeffsMax)
    {
        if (nCoeffsMin != 0)
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMin; coefLabel++)
            {
                electronicList_[coefLabel] = dem.electronicList_[coefLabel];
            }
            for (register label coefLabel=nCoeffsMin; coefLabel<nCoeffsMax; coefLabel++)
            {
                electronicList_.append(0.0);
            }
        }
        else
        {
            for (register label coefLabel=0; coefLabel<nCoeffsMax; coefLabel++)
            {
                electronicList_.append(dem.electronicList_[coefLabel]);
            }
        }
    }
    // END NEW VINCENT 11/02/2016 *********************************************
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningIn
        (
            "decoupledEnergyModesThermo<EquationOfState>::limit(const scalar T) const"
        )   << "attempt to use decoupledEnergyModesThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return min(max(T, Tlow_), Thigh_);
    }
    else
    {
        return T;
    }
}


template<class EquationOfState>
inline Foam::scalar
Foam::decoupledEnergyModesThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar
Foam::decoupledEnergyModesThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::cvCoeffsArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::decoupledCvCoeffs() const
{
    return decoupledCvCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::dynamicArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::vibrationalList() const // NEW VINCENT 24/02/2016
{
    return vibrationalList_;
}


template<class EquationOfState>
inline const typename Foam::decoupledEnergyModesThermo<EquationOfState>::dynamicArray&
Foam::decoupledEnergyModesThermo<EquationOfState>::electronicList() const
{
    return electronicList_;
}


/* Meaning of the table elements
a[0] :: cv translational
a[1] :: cv rotational
a[2] :: cv vibrational
a[3] :: cv electronic
a[4] :: cv electron
a[5] :: chemical enthalpy
a[5] :: chemical entropy
*/


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_trans
(
    const scalar p,
    const scalar Tt
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);
    return this->RR*a[0];
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_r
(
    const scalar p,
    const scalar Tt
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);
    return this->RR*a[1];
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_t
(
    const scalar p,
    const scalar Tt
) const
{
    return this->cv_trans(p, Tt) + this->cv_r(p, Tt);
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_v
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& v = vibrationalCoeffs(Tv);
    scalar tcv_v = 0.0;

    if (Tv != 0.0 and v[0] != 0.0)
    {
        if(a[2] < 0.0)
        {
            return this->RR*abs(a[2]);
        }
        else
        {
            for
            (
                label vibMode=0;
                vibMode < v.size();
                vibMode+=2
            )
            {
                scalar reducedT = v[vibMode+1] / Tv;
                tcv_v +=
                        v[vibMode] * sqr(reducedT) * exp(reducedT)
                         / sqr(exp(reducedT) - 1.0);
            }
        }
    }
    return this->RR*a[2]*tcv_v;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_v_mode
(
    const label vibMode,
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& v = vibrationalCoeffs(Tv);
    scalar tcv_v = 0.0;

    if (Tv != 0.0 and v[0] != 0.0)
    {
        if(a[2] < 0.0)
        {
            return this->RR*abs(a[2]);
        }
        else
        {
            scalar reducedT = v[2*vibMode+1] / Tv;
            tcv_v = v[2*vibMode] * sqr(reducedT) * exp(reducedT)
                     / sqr(exp(reducedT) - 1.0);
        }
    }
    return this->RR*a[2]*tcv_v;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_el
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& b = electronicCoeffs(Tv);
    scalar num1Cvel = 0.0;
    scalar num2Cvel = 0.0;
    scalar num3Cvel = 0.0;
    scalar denCvel = 0.0;

    const scalar Tel = Tv;
    if (Tel != 0.0 and a[3] != 0.0)
    {
        for
        (
            label coefLabel=2;
            coefLabel < b.size();
            coefLabel+=2
        )
        {
            num1Cvel += b[coefLabel]*sqr(b[coefLabel+1]/Tel)*exp(-b[coefLabel+1]/Tel);
            num2Cvel += b[coefLabel]*b[coefLabel+1]*exp(-b[coefLabel+1]/Tel);
            num3Cvel += b[coefLabel]*b[coefLabel+1]/sqr(Tel)*exp(-b[coefLabel+1]/Tel);
        }

        for
        (
            label coefLabel=0;
            coefLabel < b.size();
            coefLabel+=2
        )
        {
            denCvel += b[coefLabel]*exp(-b[coefLabel+1]/Tel);
        }

        return this->RR*
        (
            a[3] * (num1Cvel/denCvel - num2Cvel*num3Cvel/sqr(denCvel))
        );
    }
    else
    {
        return this->RR*0.0;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::cv_e
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);

    return this->RR*a[4];
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar Tt,
    const scalar Tv
) const
{
    return ets(p, Tt) + evs(p, Tv) + eels(p, Tv) + ees(p, Tv) + hc() + this->RR*Tt;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::ets
(
    const scalar p,
    const scalar Tt
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);
    return this->RR*a[0]*Tt + this->ers(p, Tt);
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::ers
(
    const scalar p,
    const scalar Tt
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);
    return this->RR*a[1]*Tt;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::evs_mode
(
    const label vibMode,
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& v = vibrationalCoeffs(Tv);
    scalar tev = 0.0;

    if (Tv != 0.0 and a[2] != 0.0)
    {
        if(a[2] == -1)
        {
            return this->RR*abs(a[2])*Tv;
        }
        else
        {
            scalar reducedT = v[2*vibMode+1]/Tv;
            tev = v[2*vibMode] * v[2*vibMode+1] / (exp(reducedT) - 1.0);
        }
    }

    return this->RR*a[2]*tev;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::evs
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& v = vibrationalCoeffs(Tv);
    scalar tev = 0.0;

    if (Tv != 0.0 and a[2] != 0.0)
    {
        if(a[2] == -1)
        {
            return this->RR*abs(a[2])*Tv;
        }
        else
        {
            for
            (
                label vibMode=0;
                vibMode < v.size();
                vibMode+=2
            )
            {
                scalar reducedT = v[vibMode+1]/Tv;
                tev += v[vibMode] * v[vibMode+1] / (exp(reducedT) - 1.0);
            }
        }
    }

    return this->RR*a[2]*tev;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::eels
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    const dynamicArray& b = electronicCoeffs(Tv);
    scalar numEel = 0.0;
    scalar denEel = 0.0;

    const scalar Tel = Tv;
    if (Tel != 0.0 and a[3] != 0.0)
    {
        for
        (
            label coefLabel=2;
            coefLabel < b.size();
            coefLabel+=2
        )
        {
            numEel += b[coefLabel]*b[coefLabel+1]*exp(-b[coefLabel+1]/Tel);
        }

        for
        (
            label coefLabel=0;
            coefLabel < b.size();
            coefLabel+=2
        )
        {
            denEel += b[coefLabel]*exp(-b[coefLabel+1]/Tel);
        }

        return this->RR*a[3]*numEel/denEel;
   }
   else
   {
      return this->RR*0.0;
   }

}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::ees
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    return this->RR*a[4]*Tv;
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::zetar
(
    const scalar p,
    const scalar Tt,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);

    return 2.0*a[1];
}


template<class EquationOfState>
inline Foam::scalar
Foam::decoupledEnergyModesThermo<EquationOfState>::zetav_mode
(
    const label vibMode,
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    if (Tv != 0.0 and a[2] != 0.0)
    {
        return this->evs_mode(vibMode, p, Tv) / (0.5*this->RR*Tv);
    }
    else
    {
        return 0.0;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::zetav
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    if (Tv != 0.0 and a[2] != 0.0)
    {
        return this->evs(p, Tv) / (0.5*this->RR*Tv);
    }
    else
    {
        return 0.0;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::zetael
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    if (Tv != 0.0 and a[3] != 0.0)
    {
        return this->eels(p, Tv) / (0.5*this->RR*Tv);
    }
    else
    {
        return 0.0;
    }
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::zetae
(
    const scalar p,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tv);
    if (Tv != 0.0 and a[4] != 0.0)
    {
        return this->ees(p, Tv) / (0.5*this->RR*Tv);
    }
    else
    {
        return 0.0;
    }
}


// DELETION VINCENT ***********************************************************
/*template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar Tt,
    const scalar Tv
) const
{
    return ha(p, Tt, Tv) - hc();
}*/
// END DELETION VINCENT *******************************************************


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::hc() const
{
    const scalar Tstd = advancedSpecie::Tstd;
    const cvCoeffsArray& a = cvCoeffs(Tstd);

    return this->RR*a[5];
}


template<class EquationOfState>
inline Foam::scalar Foam::decoupledEnergyModesThermo<EquationOfState>::s
(
    const scalar p,
    const scalar Tt,
    const scalar Tv
) const
{
    const cvCoeffsArray& a = cvCoeffs(Tt);
    return this->RR*a[6];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::decoupledEnergyModesThermo<EquationOfState>::operator+=
(
    const decoupledEnergyModesThermo<EquationOfState>& dem
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(dem);

    molr1 /= this->nMoles();
    scalar molr2 = dem.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, dem.Tlow_);
    Thigh_ = min(Thigh_, dem.Thigh_);

    for
    (
        register label coefLabel=0;
        coefLabel < decoupledEnergyModesThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        decoupledCvCoeffs_[coefLabel] =
            molr1*decoupledCvCoeffs_[coefLabel]
          + molr2*dem.decoupledCvCoeffs_[coefLabel];
    }

    // NEW VINCENT 23/02/2016 *************************************************
    label nCoeffsMin = min(vibrationalList_.size(), dem.vibrationalList_.size());
    label nCoeffsMax = max(vibrationalList_.size(), dem.vibrationalList_.size());

    if (nCoeffsMin != nCoeffsMax)
    {
        if (vibrationalList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList_.append(
                    molr2*dem.vibrationalList_[coefLabel]); // NEW VINCENT 23/02/2016
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList_[coefLabel] =
                    molr1*vibrationalList_[coefLabel];
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        vibrationalList_[coefLabel] =
            molr1*vibrationalList_[coefLabel]
          + molr2*dem.vibrationalList_[coefLabel];
    }
    // END NEW VINCENT 23/02/2016 *********************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(electronicList_.size(), dem.electronicList_.size());
    nCoeffsMax = max(electronicList_.size(), dem.electronicList_.size());
    /*Info << "+=/ "<< "nCoeffsMin: "<< nCoeffsMin<< ", " << "nCoeffsMax: "<< nCoeffsMax << endl;
    Info << "+=/ "<< "electronicList_.size(): "<< electronicList_.size()<< ", " << "dem.electronicList_.size(): "<< dem.electronicList_.size() << endl;*/

    if (nCoeffsMin != nCoeffsMax)
    {
        if (electronicList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList_.append(
                    molr2*dem.electronicList_[coefLabel]); // NEW VINCENT 23/02/2016
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList_[coefLabel] =
                    molr1*electronicList_[coefLabel];
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        electronicList_[coefLabel] =
            molr1*electronicList_[coefLabel]
          + molr2*dem.electronicList_[coefLabel];
    }
    // END NEW VINCENT 11/02/2016 *********************************************
}


template<class EquationOfState>
inline void Foam::decoupledEnergyModesThermo<EquationOfState>::operator-=
(
    const decoupledEnergyModesThermo<EquationOfState>& dem
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(dem);

    molr1 /= this->nMoles();
    scalar molr2 = dem.nMoles()/this->nMoles();

    Tlow_ = max(Tlow_, dem.Tlow_);
    Thigh_ = min(Thigh_, dem.Thigh_);

    for
    (
        register label coefLabel=0;
        coefLabel<decoupledEnergyModesThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        decoupledCvCoeffs_[coefLabel] =
            molr1*decoupledCvCoeffs_[coefLabel]
          - molr2*dem.decoupledCvCoeffs_[coefLabel];
    }

    // NEW VINCENT 23/02/2016 *************************************************
    label nCoeffsMin = min(vibrationalList_.size(), dem.vibrationalList_.size());
    label nCoeffsMax = max(vibrationalList_.size(), dem.vibrationalList_.size());

    if (nCoeffsMin != nCoeffsMax)
    {
        if (vibrationalList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList_.append(
                    - molr2*dem.vibrationalList_[coefLabel]); // NEW VINCENT 23/02/2016
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList_[coefLabel] =
                    molr1*vibrationalList_[coefLabel];
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        vibrationalList_[coefLabel] =
            molr1*vibrationalList_[coefLabel]
          - molr2*dem.vibrationalList_[coefLabel];
    }
    // NEW VINCENT 23/02/2016 *************************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(electronicList_.size(), dem.electronicList_.size());
    nCoeffsMax = max(electronicList_.size(), dem.electronicList_.size());
    //Info << "-=/ "<< "nCoeffsMin: "<< nCoeffsMin<< ", " << "nCoeffsMax: "<< nCoeffsMax << endl;
    //Info << "-=/ "<< "electronicList_.size(): "<< electronicList_.size()<< ", " << "dem.electronicList_.size(): "<< dem.electronicList_.size() << endl;

    if (nCoeffsMin != nCoeffsMax)
    {
        if (electronicList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList_.append(
                    - molr2*dem.electronicList_[coefLabel]); // NEW VINCENT 23/02/2016
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList_[coefLabel] =
                    molr1*electronicList_[coefLabel];
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        electronicList_[coefLabel] =
            molr1*electronicList_[coefLabel]
          - molr2*dem.electronicList_[coefLabel];
    }
    // NEW VINCENT 11/02/2016 *************************************************
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::decoupledEnergyModesThermo<EquationOfState> Foam::operator+
(
    const decoupledEnergyModesThermo<EquationOfState>& dem1,
    const decoupledEnergyModesThermo<EquationOfState>& dem2
)
{
    EquationOfState eofs = dem1;
    eofs += dem2;

    scalar molr1 = dem1.nMoles()/eofs.nMoles();
    scalar molr2 = dem2.nMoles()/eofs.nMoles();

    typename decoupledEnergyModesThermo<EquationOfState>::cvCoeffsArray decoupledCvCoeffs;
    typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray vibrationalList; // NEW VINCENT 24/02/2016
    typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray electronicList;

    for
    (
        register label coefLabel=0;
        coefLabel<decoupledEnergyModesThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        decoupledCvCoeffs[coefLabel] =
            molr1*dem1.decoupledCvCoeffs_[coefLabel]
          + molr2*dem2.decoupledCvCoeffs_[coefLabel];
    }

    // NEW VINCENT 24/02/2016 *************************************************
    label nCoeffsMin = min(dem1.vibrationalList_.size(), dem2.vibrationalList_.size());
    label nCoeffsMax = max(dem1.vibrationalList_.size(), dem2.vibrationalList_.size());

    if (nCoeffsMin != nCoeffsMax)
    {
        if (dem1.vibrationalList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList.append(
                    molr2*dem2.vibrationalList_[coefLabel]);
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList.append(
                    molr1*dem1.vibrationalList_[coefLabel]);
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        vibrationalList.append(
            molr1*dem1.vibrationalList_[coefLabel]
          + molr2*dem2.vibrationalList_[coefLabel]);
    }
    // NEW VINCENT 24/02/2016 *************************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(dem1.electronicList_.size(), dem2.electronicList_.size());
    nCoeffsMax = max(dem1.electronicList_.size(), dem2.electronicList_.size());
    //Info << "+/ "<< "nCoeffsMin: "<< nCoeffsMin<< ", " << "nCoeffsMax: "<< nCoeffsMax << endl;
    //Info << "+/ "<< "dem1.electronicList_.size(): "<< dem1.electronicList_.size()<< ", " << "dem2.electronicList_.size(): "<< dem2.electronicList_.size() << endl;

    if (nCoeffsMin != nCoeffsMax)
    {
        if (dem1.electronicList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList.append(
                    molr2*dem2.electronicList_[coefLabel]);
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList.append(
                    molr1*dem1.electronicList_[coefLabel]);
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        electronicList.append(
            molr1*dem1.electronicList_[coefLabel]
          + molr2*dem2.electronicList_[coefLabel]);
    }
    // NEW VINCENT 11/02/2016 *************************************************

    return decoupledEnergyModesThermo<EquationOfState>
    (
        eofs,
        max(dem1.Tlow_, dem2.Tlow_),
        min(dem1.Thigh_, dem2.Thigh_),
        decoupledCvCoeffs,
        vibrationalList, // NEW VINCENT 24/02/2016
        electronicList
    );
}


template<class EquationOfState>
inline Foam::decoupledEnergyModesThermo<EquationOfState> Foam::operator-
(
    const decoupledEnergyModesThermo<EquationOfState>& dem1,
    const decoupledEnergyModesThermo<EquationOfState>& dem2
)
{
    EquationOfState eofs = dem1;
    eofs -= dem2;

    scalar molr1 = dem1.nMoles()/eofs.nMoles();
    scalar molr2 = dem2.nMoles()/eofs.nMoles();

    typename decoupledEnergyModesThermo<EquationOfState>::cvCoeffsArray decoupledCvCoeffs;
    typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray vibrationalList; // NEW VINCENT 24/02/2016
    typename decoupledEnergyModesThermo<EquationOfState>::dynamicArray electronicList;

    for
    (
        register label coefLabel=0;
        coefLabel<decoupledEnergyModesThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        decoupledCvCoeffs[coefLabel] =
            molr1*dem1.decoupledCvCoeffs_[coefLabel]
          - molr2*dem2.decoupledCvCoeffs_[coefLabel];
    }

    // NEW VINCENT 24/02/2016 *************************************************
    label nCoeffsMin = min(dem1.vibrationalList_.size(), dem2.vibrationalList_.size());
    label nCoeffsMax = max(dem1.vibrationalList_.size(), dem2.vibrationalList_.size());

    if (nCoeffsMin != nCoeffsMax)
    {
        if (dem1.vibrationalList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList.append(
                    - molr2*dem2.vibrationalList_[coefLabel]);
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                vibrationalList.append(
                    molr1*dem1.vibrationalList_[coefLabel]);
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        vibrationalList.append(
            molr1*dem1.vibrationalList_[coefLabel]
          - molr2*dem2.vibrationalList_[coefLabel]);
    }
    // NEW VINCENT 24/02/2016 *************************************************

    // NEW VINCENT 11/02/2016 *************************************************
    nCoeffsMin = min(dem1.electronicList_.size(), dem2.electronicList_.size());
    nCoeffsMax = max(dem1.electronicList_.size(), dem2.electronicList_.size());
    //Info << "-/ "<< "nCoeffsMin: "<< nCoeffsMin<< ", " << "nCoeffsMax: "<< nCoeffsMax << endl;
    //Info << "-/ "<< "dem1.electronicList_.size(): "<< dem1.electronicList_.size()<< ", " << "dem2.electronicList_.size(): "<< dem2.electronicList_.size() << endl;

    if (nCoeffsMin != nCoeffsMax)
    {
        if (dem1.electronicList_.size() < nCoeffsMax)
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList.append(
                    - molr2*dem2.electronicList_[coefLabel]);
            }
        }
        else
        {
            for
            (
                register label coefLabel=nCoeffsMin;
                coefLabel < nCoeffsMax;
                coefLabel++
            )
            {
                electronicList.append(
                    molr1*dem1.electronicList_[coefLabel]);
            }
        }
    }

    for
    (
        register label coefLabel=0;
        coefLabel < nCoeffsMin;
        coefLabel++
    )
    {
        electronicList.append(
            molr1*dem1.electronicList_[coefLabel]
          - molr2*dem2.electronicList_[coefLabel]);
    }
    // NEW VINCENT 11/02/2016 *************************************************

    return decoupledEnergyModesThermo<EquationOfState>
    (
        eofs,
        max(dem1.Tlow_, dem2.Tlow_),
        min(dem1.Thigh_, dem2.Thigh_),
        decoupledCvCoeffs,
        vibrationalList, // NEW VINCENT 24/02/2016
        electronicList
    );
}


template<class EquationOfState>
inline Foam::decoupledEnergyModesThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const decoupledEnergyModesThermo<EquationOfState>& dem
)
{
    return decoupledEnergyModesThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(dem),
        dem.Tlow_,
        dem.Thigh_,
        dem.decoupledCvCoeffs_,
        dem.vibrationalList_, // NEW VINCENT 24/02/2016
        dem.electronicList_
    );
}


template<class EquationOfState>
inline Foam::decoupledEnergyModesThermo<EquationOfState> Foam::operator==
(
    const decoupledEnergyModesThermo<EquationOfState>& dem1,
    const decoupledEnergyModesThermo<EquationOfState>& dem2
)
{
    return dem2 - dem1;
}


// ************************************************************************* //
