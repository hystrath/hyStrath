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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constantTransport<Thermo>::constantTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar kappa,
    const scalar kappave
)
:
    Thermo(t),
    mu_(mu),
    kappa_(kappa),
    kappave_(kappave)
{}


template<class Thermo>
inline Foam::constantTransport<Thermo>::constantTransport
(
    const word& name,
    const constantTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    kappa_(ct.kappa_),
    kappave_(ct.kappave_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constantTransport<Thermo> >
Foam::constantTransport<Thermo>::clone() const
{
    return autoPtr<constantTransport<Thermo> >
    (
        new constantTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constantTransport<Thermo> >
Foam::constantTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<constantTransport<Thermo> >
    (
        new constantTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::constantTransport<Thermo> >
Foam::constantTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<constantTransport<Thermo> >
    (
        new constantTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::kappatrans
(
    const scalar p, const scalar Tt
) const
{
    return kappa_;
}


template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::kappar
(
    const scalar p, const scalar Tt
) const
{
    return 0.0;
}


template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::kappatr
(
    const scalar p, const scalar Tt
) const
{
    return kappatrans(p, Tt) + kappar(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::kappave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    return kappave_;
}


template<class Thermo>
inline Foam::scalar Foam::constantTransport<Thermo>::alphatr
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
inline Foam::scalar Foam::constantTransport<Thermo>::alphave
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
inline Foam::constantTransport<Thermo>&
Foam::constantTransport<Thermo>::operator=
(
    const constantTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    kappa_ = ct.kappa_;
    kappave_ = ct.kappave_;

    return *this;
}


template<class Thermo>
inline void Foam::constantTransport<Thermo>::operator+=
(
    const constantTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ + molr2*st.mu_;
    kappa_ = molr1*kappa_ + molr2*st.kappa_;
    kappave_ = molr1*kappave_ + molr2*st.kappave_;
}


template<class Thermo>
inline void Foam::constantTransport<Thermo>::operator-=
(
    const constantTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ - molr2*st.mu_;
    kappa_ = molr1*kappa_ - molr2*st.kappa_;
    kappave_ = molr1*kappave_ - molr2*st.kappave_;
}


template<class Thermo>
inline void Foam::constantTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::constantTransport<Thermo> Foam::operator+
(
    const constantTransport<Thermo>& ct1,
    const constantTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constantTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ + molr2*ct2.mu_,
        molr1*ct1.kappa_ + molr2*ct2.kappa_,
        molr1*ct1.kappave_ + molr2*ct2.kappave_
    );
}


template<class Thermo>
inline Foam::constantTransport<Thermo> Foam::operator-
(
    const constantTransport<Thermo>& ct1,
    const constantTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constantTransport<Thermo>
    (
        t,
        molr1*ct1.mu_ - molr2*ct2.mu_,
        molr1*ct1.kappa_ - molr2*ct2.kappa_,
        molr1*ct1.kappave_ - molr2*ct2.kappave_
    );
}


template<class Thermo>
inline Foam::constantTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constantTransport<Thermo>& ct
)
{
    return constantTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.kappa_,
        ct.kappave_
    );
}


template<class Thermo>
inline Foam::constantTransport<Thermo> Foam::operator==
(
    const constantTransport<Thermo>& ct1,
    const constantTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
