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

#include "fundamentalConstants.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport
(
    const Thermo& t,
    const scalar dref,
    const scalar omega,
    const scalar eta_s
)
:
    Thermo(t),
    dref_(dref),
    omega_(omega),
    eta_s_(eta_s)
{}


template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo>::powerLawEuckenTransport
(
    const word& name,
    const powerLawEuckenTransport& pet
)
:
    Thermo(name, pet),
    dref_(pet.dref_),
    omega_(pet.omega_),
    eta_s_(pet.eta_s_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawEuckenTransport<Thermo> >
Foam::powerLawEuckenTransport<Thermo>::clone() const
{
    return autoPtr<powerLawEuckenTransport<Thermo> >
    (
        new powerLawEuckenTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawEuckenTransport<Thermo> >
Foam::powerLawEuckenTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<powerLawEuckenTransport<Thermo> >
    (
        new powerLawEuckenTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawEuckenTransport<Thermo> >
Foam::powerLawEuckenTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<powerLawEuckenTransport<Thermo> >
    (
        new powerLawEuckenTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::mu
(
    const scalar p,
    const scalar Tt
) const
{
    scalar Tref = 273.0;
    scalar mu_ref = 15.0
        *sqrt(constant::mathematical::pi*constant::physicoChemical::k.value()
        *(this->W()*1.0e-3/constant::physicoChemical::NA.value())*Tref)
        /(constant::mathematical::twoPi*pow(dref_, 2.0)*(5.0-2.0*omega_)
            *(7.0-2.0*omega_));
    
    return mu_ref*pow(Tt/Tref, omega_);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::kappatrans
(
    const scalar p, const scalar Tt
) const
{
    return 5.0/2.0*mu(p, Tt)*this->Cv_trans(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::kappar
(
    const scalar p, const scalar Tt
) const
{
    return mu(p, Tt)*this->Cv_r(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::kappatr
(
    const scalar p, const scalar Tt
) const
{
    return kappatrans(p, Tt) + kappar(p, Tt);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::kappave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    return eta_s_*mu(p, Tt)*this->Cv_vel(p, Tve);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::alphatr
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
inline Foam::scalar Foam::powerLawEuckenTransport<Thermo>::alphave
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
inline Foam::powerLawEuckenTransport<Thermo>&
Foam::powerLawEuckenTransport<Thermo>::operator=
(
    const powerLawEuckenTransport<Thermo>& pet
)
{
    Thermo::operator=(pet);

    dref_ = pet.dref_;
    omega_ = pet.omega_;
    eta_s_ = pet.eta_s_;

    return *this;
}


template<class Thermo>
inline void Foam::powerLawEuckenTransport<Thermo>::operator+=
(
    const powerLawEuckenTransport<Thermo>& pet
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(pet);

    molr1 /= this->nMoles();
    scalar molr2 = pet.nMoles()/this->nMoles();

    dref_ = molr1*dref_ + molr2*pet.dref_;
    omega_ = molr1*omega_ + molr2*pet.omega_;
    eta_s_ = molr1*eta_s_ + molr2*pet.eta_s_;
}


template<class Thermo>
inline void Foam::powerLawEuckenTransport<Thermo>::operator-=
(
    const powerLawEuckenTransport<Thermo>& pet
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(pet);

    molr1 /= this->nMoles();
    scalar molr2 = pet.nMoles()/this->nMoles();

    dref_ = molr1*dref_ - molr2*pet.dref_;
    omega_ = molr1*omega_ - molr2*pet.omega_;
    eta_s_ = molr1*eta_s_ - molr2*pet.eta_s_;
}


template<class Thermo>
inline void Foam::powerLawEuckenTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo> Foam::operator+
(
    const powerLawEuckenTransport<Thermo>& pet1,
    const powerLawEuckenTransport<Thermo>& pet2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pet1) + static_cast<const Thermo&>(pet2)
    );

    scalar molr1 = pet1.nMoles()/t.nMoles();
    scalar molr2 = pet2.nMoles()/t.nMoles();

    return powerLawEuckenTransport<Thermo>
    (
        t,
        molr1*pet1.dref_ + molr2*pet2.dref_,
        molr1*pet1.omega_ + molr2*pet2.omega_,
        molr1*pet1.eta_s_ + molr2*pet2.eta_s_
    );
}


template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo> Foam::operator-
(
    const powerLawEuckenTransport<Thermo>& pet1,
    const powerLawEuckenTransport<Thermo>& pet2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pet1) - static_cast<const Thermo&>(pet2)
    );

    scalar molr1 = pet1.nMoles()/t.nMoles();
    scalar molr2 = pet2.nMoles()/t.nMoles();

    return powerLawEuckenTransport<Thermo>
    (
        t,
        molr1*pet1.dref_ - molr2*pet2.dref_,
        molr1*pet1.omega_ - molr2*pet2.omega_,
        molr1*pet1.eta_s_ - molr2*pet2.eta_s_
    );
}


template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo> Foam::operator*
(
    const scalar s,
    const powerLawEuckenTransport<Thermo>& pet
)
{
    return powerLawEuckenTransport<Thermo>
    (
        s*static_cast<const Thermo&>(pet),
        pet.dref_,
        pet.omega_,
        pet.eta_s_
    );
}


template<class Thermo>
inline Foam::powerLawEuckenTransport<Thermo> Foam::operator==
(
    const powerLawEuckenTransport<Thermo>& pet1,
    const powerLawEuckenTransport<Thermo>& pet2
)
{
    return pet2 - pet1;
}


// ************************************************************************* //
