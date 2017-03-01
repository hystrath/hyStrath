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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport
(
    const Thermo& t,
    const scalar Ak,
    const scalar Bk,
    const scalar Ck
)
:
    Thermo(t),
    Ak_(Ak),
    Bk_(Bk),
    Ck_(Ck)
{}


template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo>::BlottnerEuckenTransport
(
    const word& name,
    const BlottnerEuckenTransport& bet
)
:
    Thermo(name, bet),
    Ak_(bet.Ak_),
    Bk_(bet.Bk_),
    Ck_(bet.Ck_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::BlottnerEuckenTransport<Thermo> >
Foam::BlottnerEuckenTransport<Thermo>::clone() const
{
    return autoPtr<BlottnerEuckenTransport<Thermo> >
    (
        new BlottnerEuckenTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::BlottnerEuckenTransport<Thermo> >
Foam::BlottnerEuckenTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<BlottnerEuckenTransport<Thermo> >
    (
        new BlottnerEuckenTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::BlottnerEuckenTransport<Thermo> >
Foam::BlottnerEuckenTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<BlottnerEuckenTransport<Thermo> >
    (
        new BlottnerEuckenTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::mu
(
    const scalar p,
    const scalar Tt
) const
{
    return 0.1*::exp((Ak_*::log(Tt) + Bk_)*::log(Tt) + Ck_);
}


template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::kappatrans
(
    const scalar p, const scalar Tt
) const
{
    return 5.0/2.0*mu(p, Tt)*this->Cv_trans(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::kappar
(
    const scalar p, const scalar Tt
) const
{
    return mu(p, Tt)*this->Cv_r(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::kappatr
(
    const scalar p, const scalar Tt
) const
{
    return kappatrans(p, Tt) + kappar(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::kappave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    const scalar eta_s = 1.2; // Hard-sphere: (7-2*omega_m)/5 with omega_m = 0.5
    return eta_s*mu(p, Tt)*this->Cv_vel(p, Tve);
}


template<class Thermo>
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::alphatr
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
inline Foam::scalar Foam::BlottnerEuckenTransport<Thermo>::alphave
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
inline Foam::BlottnerEuckenTransport<Thermo>&
Foam::BlottnerEuckenTransport<Thermo>::operator=
(
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    Thermo::operator=(bet);

    Ak_ = bet.Ak_;
    Bk_ = bet.Bk_;
    Ck_ = bet.Ck_;

    return *this;
}


template<class Thermo>
inline void Foam::BlottnerEuckenTransport<Thermo>::operator+=
(
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(bet);

    molr1 /= this->nMoles();
    scalar molr2 = bet.nMoles()/this->nMoles();

    Ak_ = molr1*Ak_ + molr2*bet.Ak_;
    Bk_ = molr1*Bk_ + molr2*bet.Bk_;
    Ck_ = molr1*Ck_ + molr2*bet.Ck_;
}


template<class Thermo>
inline void Foam::BlottnerEuckenTransport<Thermo>::operator-=
(
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(bet);

    molr1 /= this->nMoles();
    scalar molr2 = bet.nMoles()/this->nMoles();

    Ak_ = molr1*Ak_ - molr2*bet.Ak_;
    Bk_ = molr1*Bk_ - molr2*bet.Bk_;
    Ck_ = molr1*Ck_ - molr2*bet.Ck_;
}


template<class Thermo>
inline void Foam::BlottnerEuckenTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo> Foam::operator+
(
    const BlottnerEuckenTransport<Thermo>& bet1,
    const BlottnerEuckenTransport<Thermo>& bet2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(bet1) + static_cast<const Thermo&>(bet2)
    );

    scalar molr1 = bet1.nMoles()/t.nMoles();
    scalar molr2 = bet2.nMoles()/t.nMoles();

    return BlottnerEuckenTransport<Thermo>
    (
        t,
        molr1*bet1.Ak_ + molr2*bet2.Ak_,
        molr1*bet1.Bk_ + molr2*bet2.Bk_,
        molr1*bet1.Ck_ + molr2*bet2.Ck_
    );
}


template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo> Foam::operator-
(
    const BlottnerEuckenTransport<Thermo>& bet1,
    const BlottnerEuckenTransport<Thermo>& bet2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(bet1) - static_cast<const Thermo&>(bet2)
    );

    scalar molr1 = bet1.nMoles()/t.nMoles();
    scalar molr2 = bet2.nMoles()/t.nMoles();

    return BlottnerEuckenTransport<Thermo>
    (
        t,
        molr1*bet1.Ak_ - molr2*bet2.Ak_,
        molr1*bet1.Bk_ - molr2*bet2.Bk_,
        molr1*bet1.Ck_ - molr2*bet2.Ck_
    );
}


template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo> Foam::operator*
(
    const scalar s,
    const BlottnerEuckenTransport<Thermo>& bet
)
{
    return BlottnerEuckenTransport<Thermo>
    (
        s*static_cast<const Thermo&>(bet),
        bet.Ak_,
        bet.Bk_,
        bet.Ck_
    );
}


template<class Thermo>
inline Foam::BlottnerEuckenTransport<Thermo> Foam::operator==
(
    const BlottnerEuckenTransport<Thermo>& bet1,
    const BlottnerEuckenTransport<Thermo>& bet2
)
{
    return bet2 - bet1;
}


// ************************************************************************* //
