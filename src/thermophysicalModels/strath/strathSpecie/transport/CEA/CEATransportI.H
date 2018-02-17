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
inline Foam::CEATransport<Thermo>::CEATransport
(
    const Thermo& t,
    const DynamicList<scalar> temp,
    const DynamicList<CEATransportArray> mu,
    const DynamicList<CEATransportArray> kappa
)
:
    Thermo(t),
    temp_(temp),
    mu_(mu),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::CEATransport<Thermo>::CEATransport
(
    const word& name,
    const CEATransport& ceat
)
:
    Thermo(name, ceat),
    temp_(ceat.temp_),
    mu_(ceat.mu_),
    kappa_(ceat.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::CEATransport<Thermo> >
Foam::CEATransport<Thermo>::clone() const
{
    return autoPtr<CEATransport<Thermo> >
    (
        new CEATransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CEATransport<Thermo> >
Foam::CEATransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<CEATransport<Thermo> >
    (
        new CEATransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CEATransport<Thermo> >
Foam::CEATransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<CEATransport<Thermo> >
    (
        new CEATransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline const typename Foam::CEATransport<Thermo>::CEATransportArray&
Foam::CEATransport<Thermo>::viscosityCoeffs
(
    const scalar T
) const
{
    for(int Ti=0; Ti < temp_.size(); Ti++)
    {
        if(T < temp_[Ti+1])
        {
            return mu_[Ti];
        }
    }
    
    return mu_[temp_.size()-1];
}


template<class Thermo>
inline const typename Foam::CEATransport<Thermo>::CEATransportArray&
Foam::CEATransport<Thermo>::thermalConductivityCoeffs
(
    const scalar T
) const
{
    for(int Ti=0; Ti < temp_.size(); Ti++)
    {
        if(T < temp_[Ti+1])
        {
            return kappa_[Ti];
        }
    }
    
    return kappa_[temp_.size()-1];
}



template<class Thermo>
inline Foam::scalar Foam::CEATransport<Thermo>::mu
(
    const scalar p,
    const scalar Tt
) const
{
    const CEATransportArray& muCoeffs = viscosityCoeffs(Tt);
    
    //NB: The coefficients in the dictionary tables are given in microPoise (cf Eilmer3, CEA-viscosity.cxx)
    return 1.0e-7*Foam::exp(muCoeffs[0]*::log(Tt) + muCoeffs[1]/Tt + muCoeffs[2]/Foam::sqr(Tt) + muCoeffs[3]);
}


template<class Thermo>
inline Foam::scalar Foam::CEATransport<Thermo>::kappatr
(
    const scalar p, const scalar Tt
) const
{
    const CEATransportArray& kappaCoeffs = thermalConductivityCoeffs(Tt);
    
    //NB: The coefficients in the dictionary tables are given in microWatts/(cm.K) (cf Eilmer3, CEA-thermal-conductivity.cxx)
    return 1.0e-4*Foam::exp(kappaCoeffs[0]*::log(Tt) + kappaCoeffs[1]/Tt + kappaCoeffs[2]/Foam::sqr(Tt) + kappaCoeffs[3]);
}


template<class Thermo>
inline Foam::scalar Foam::CEATransport<Thermo>::kappave
(
    const scalar p, const scalar Tt, const scalar Tve
) const
{
    const scalar eta_s = 1.2;
    return eta_s*mu(p, Tt)*this->Cv_vel(p, Tve);
}


template<class Thermo>
inline Foam::scalar Foam::CEATransport<Thermo>::alphatr
(
    const scalar p, const scalar Tt
) const
{
    if (this->Cv_t(p, Tt) != 0.0)
    {
        //Info << "Hi 2.0" << endl;
        //Info << "kappatr(p, Tt): " << kappatr(p, Tt) << endl;
        //Info << "this->Cv_t(p, Tt): " << this->Cv_t(p, Tt) << endl;
        return kappatr(p, Tt)/this->Cv_t(p, Tt);
    }
    else
    {
        //Info << "Hi 3.0" << endl;
        return 0.0;
    } 
}


template<class Thermo>
inline Foam::scalar Foam::CEATransport<Thermo>::alphave
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
inline Foam::CEATransport<Thermo>&
Foam::CEATransport<Thermo>::operator=
(
    const CEATransport<Thermo>& ceat
)
{
    Thermo::operator=(ceat);

    temp_ = ceat.temp_;
    mu_ = ceat.mu_;
    kappa_ = ceat.kappa_;

    return *this;
}


template<class Thermo>
inline void Foam::CEATransport<Thermo>::operator+=
(
    const CEATransport<Thermo>& ceat
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(ceat);

    molr1 /= this->nMoles();
    scalar molr2 = ceat.nMoles()/this->nMoles();

    temp_ = molr1*temp_ + molr2*ceat.temp_;

    for
    (
        register label coefLabel=0;
        coefLabel<ceat.mu_.size();
        coefLabel++
    )
    {
        for
        (
            register label coefLabelj=0;
            coefLabelj<CEATransport<Thermo>::nCoeffs_;
            coefLabelj++
        )
        {
            mu_[coefLabel][coefLabelj] =
                molr1*mu_[coefLabel][coefLabelj]
              + molr2*ceat.mu_[coefLabel][coefLabelj];

            kappa_[coefLabel][coefLabelj] =
                molr1*kappa_[coefLabel][coefLabelj]
              + molr2*ceat.kappa_[coefLabel][coefLabelj];
        }
    }
    
    
    //mu_ = molr1*mu_ + molr2*ceat.mu_;
    //kappa_ = molr1*kappa_ + molr2*ceat.kappa_;
}


template<class Thermo>
inline void Foam::CEATransport<Thermo>::operator-=
(
    const CEATransport<Thermo>& ceat
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(ceat);

    molr1 /= this->nMoles();
    scalar molr2 = ceat.nMoles()/this->nMoles();

    temp_ = molr1*temp_ - molr2*ceat.temp_;

    for
    (
        register label coefLabel=0;
        coefLabel<ceat.mu_.size();
        coefLabel++
    )
    {
        for
        (
            register label coefLabelj=0;
            coefLabelj<CEATransport<Thermo>::nCoeffs_;
            coefLabelj++
        )
        {
            mu_[coefLabel][coefLabelj] =
                molr1*mu_[coefLabel][coefLabelj]
              - molr2*ceat.mu_[coefLabel][coefLabelj];

            kappa_[coefLabel][coefLabelj] =
                molr1*kappa_[coefLabel][coefLabelj]
              - molr2*ceat.kappa_[coefLabel][coefLabelj];
        }
    }
    //mu_ = molr1*mu_ - molr2*ceat.mu_;
    //kappa_ = molr1*kappa_ - molr2*ceat.kappa_;
}


template<class Thermo>
inline void Foam::CEATransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CEATransport<Thermo> Foam::operator+
(
    const CEATransport<Thermo>& ceat1,
    const CEATransport<Thermo>& ceat2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ceat1) + static_cast<const Thermo&>(ceat2)
    );

    scalar molr1 = ceat1.nMoles()/t.nMoles();
    scalar molr2 = ceat2.nMoles()/t.nMoles();

    DynamicList<FixedList<scalar, Foam::CEATransport<Thermo>::nCoeffs_> > mu, kappa;
    /*for
    (
        register label coefLabel=0;
        coefLabel<ceat1.mu_.size();
        coefLabel++
    )
    {
        for
        (
            register label coefLabelj=0;
            coefLabelj<CEATransport<Thermo>::nCoeffs_;
            coefLabelj++
        )
        {
            mu[coefLabel][coefLabelj] =
                molr1*ceat1.mu_[coefLabel][coefLabelj]
              + molr2*ceat2.mu_[coefLabel][coefLabelj];

            kappa[coefLabel][coefLabelj] =
                molr1*ceat1.kappa_[coefLabel][coefLabelj]
              + molr2*ceat2.kappa_[coefLabel][coefLabelj];
        }
    }*/
    
    return CEATransport<Thermo>
    (
        t,
        molr1*ceat1.temp_ + molr2*ceat2.temp_,
        mu,
        kappa
    );
}


template<class Thermo>
inline Foam::CEATransport<Thermo> Foam::operator-
(
    const CEATransport<Thermo>& ceat1,
    const CEATransport<Thermo>& ceat2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ceat1) - static_cast<const Thermo&>(ceat2)
    );

    scalar molr1 = ceat1.nMoles()/t.nMoles();
    scalar molr2 = ceat2.nMoles()/t.nMoles();
    
    DynamicList<FixedList<scalar, Foam::CEATransport<Thermo>::nCoeffs_> > mu, kappa;
    /*for
    (
        register label coefLabel=0;
        coefLabel<ceat1.mu_.size();
        coefLabel++
    )
    {
        for
        (
            register label coefLabelj=0;
            coefLabelj<CEATransport<Thermo>::nCoeffs_;
            coefLabelj++
        )
        {
            mu[coefLabel][coefLabelj] =
                molr1*ceat1.mu_[coefLabel][coefLabelj]
              - molr2*ceat2.mu_[coefLabel][coefLabelj];

            kappa[coefLabel][coefLabelj] =
                molr1*ceat1.kappa_[coefLabel][coefLabelj]
              - molr2*ceat2.kappa_[coefLabel][coefLabelj];
        }
    }*/

    return CEATransport<Thermo>
    (
        t,
        molr1*ceat1.temp_ - molr2*ceat2.temp_,
        mu,
        kappa
    );
}


template<class Thermo>
inline Foam::CEATransport<Thermo> Foam::operator*
(
    const scalar s,
    const CEATransport<Thermo>& ceat
)
{
    return CEATransport<Thermo>
    (
        s*static_cast<const Thermo&>(ceat),
        ceat.temp_,
        ceat.mu_,
        ceat.kappa_
    );
}


template<class Thermo>
inline Foam::CEATransport<Thermo> Foam::operator==
(
    const CEATransport<Thermo>& ceat1,
    const CEATransport<Thermo>& ceat2
)
{
    return ceat2 - ceat1;
}


// ************************************************************************* //
