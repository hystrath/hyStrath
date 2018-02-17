/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "advancedSpecie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
inline void Foam::SutherlandEuckenTransport<Thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts)
{}


template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    Thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo>::SutherlandEuckenTransport
(
    const word& name,
    const SutherlandEuckenTransport& st
)
:
    Thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::SutherlandEuckenTransport<Thermo> >
Foam::SutherlandEuckenTransport<Thermo>::clone() const
{
    return autoPtr<SutherlandEuckenTransport<Thermo> >
    (
        new SutherlandEuckenTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::SutherlandEuckenTransport<Thermo> >
Foam::SutherlandEuckenTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<SutherlandEuckenTransport<Thermo> >
    (
        new SutherlandEuckenTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::SutherlandEuckenTransport<Thermo> >
Foam::SutherlandEuckenTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<SutherlandEuckenTransport<Thermo> >
    (
        new SutherlandEuckenTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::kappatrans
(
    const scalar p, const scalar Tt
) const
{
    return 5.0/2.0*mu(p, Tt)*this->Cv_trans(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::kappar
(
    const scalar p, const scalar Tt
) const
{
    return mu(p, Tt)*this->Cv_r(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::kappatr
(
    const scalar p, const scalar Tt
) const
{
    scalar Cv_ = this->Cv_t(p, Tt);
    return mu(p, Tt)*Cv_*(1.32 + 1.77*this->R()/Cv_);
    //return kappatrans(p, Tt) + kappar(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::kappave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    const scalar eta_s = 1.2;
    return eta_s*mu(p, Tt)*this->Cv_vel(p, Tve);
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::alphatr
(
    const scalar p, const scalar Tt
) const
{
    if (this->Cv_t(p, Tt) != 0.0)
    {
        return kappatr(p, Tt)/this->Cv_t(p, Tt);
    }
    else
    {
        return 0.0;
    }    
}


template<class Thermo>
inline Foam::scalar Foam::SutherlandEuckenTransport<Thermo>::alphave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    if (this->Cv_vel(p, Tve) != 0.0)
    {
        return kappave(p, Tt, Tve)/this->Cv_vel(p, Tve);
    }
    else
    {
        return 0.0;
    }    
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo>&
Foam::SutherlandEuckenTransport<Thermo>::operator=
(
    const SutherlandEuckenTransport<Thermo>& st
)
{
    Thermo::operator=(st);

    As_ = st.As_;
    Ts_ = st.Ts_;

    return *this;
}


template<class Thermo>
inline void Foam::SutherlandEuckenTransport<Thermo>::operator+=
(
    const SutherlandEuckenTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    As_ = molr1*As_ + molr2*st.As_;
    Ts_ = molr1*Ts_ + molr2*st.Ts_;
}


template<class Thermo>
inline void Foam::SutherlandEuckenTransport<Thermo>::operator-=
(
    const SutherlandEuckenTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    As_ = molr1*As_ - molr2*st.As_;
    Ts_ = molr1*Ts_ - molr2*st.Ts_;
}


template<class Thermo>
inline void Foam::SutherlandEuckenTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo> Foam::operator+
(
    const SutherlandEuckenTransport<Thermo>& st1,
    const SutherlandEuckenTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return SutherlandEuckenTransport<Thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo> Foam::operator-
(
    const SutherlandEuckenTransport<Thermo>& st1,
    const SutherlandEuckenTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return SutherlandEuckenTransport<Thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_
    );
}


template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo> Foam::operator*
(
    const scalar s,
    const SutherlandEuckenTransport<Thermo>& st
)
{
    return SutherlandEuckenTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.As_,
        st.Ts_
    );
}


template<class Thermo>
inline Foam::SutherlandEuckenTransport<Thermo> Foam::operator==
(
    const SutherlandEuckenTransport<Thermo>& st1,
    const SutherlandEuckenTransport<Thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
