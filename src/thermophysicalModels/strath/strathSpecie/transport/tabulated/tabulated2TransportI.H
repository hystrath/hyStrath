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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulated2Transport<Thermo>::tabulated2Transport
(
    const Thermo& t,
    const interpolation2DTable<scalar> mu,
    const interpolation2DTable<scalar> kappa
)
:
    Thermo(t),
    mu_(mu),
    kappa_(kappa)
{
    //mu_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
    //kappa_.outOfBounds(interpolation2DTable<scalar>::EXTRAPOLATE);
}


template<class Thermo>
inline Foam::tabulated2Transport<Thermo>::tabulated2Transport
(
    const word& name,
    const tabulated2Transport& tt
)
:
    Thermo(name, tt),
    mu_(tt.mu_),
    kappa_(tt.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2Transport<Thermo> >
Foam::tabulated2Transport<Thermo>::clone() const
{
    return autoPtr<tabulated2Transport<Thermo> >
    (
        new tabulated2Transport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2Transport<Thermo> >
Foam::tabulated2Transport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<tabulated2Transport<Thermo> >
    (
        new tabulated2Transport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2Transport<Thermo> >
Foam::tabulated2Transport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulated2Transport<Thermo> >
    (
        new tabulated2Transport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulated2Transport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_(T, p);
}


template<class Thermo>
inline Foam::scalar Foam::tabulated2Transport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}
/*template<class Thermo>
inline Foam::scalar Foam::tabulated2Transport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_(T, p);
}*/


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulated2Transport<Thermo>& Foam::tabulated2Transport<Thermo>::operator=
(
    const tabulated2Transport<Thermo>& tt
)
{
    Thermo::operator=(tt);

    mu_ = tt.mu_;
    kappa_ = tt.kappa_;

    return *this;
}


template<class Thermo>
inline void Foam::tabulated2Transport<Thermo>::operator+=
(
    const tabulated2Transport<Thermo>& tt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(tt);
    /*
    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ + molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    */
}


template<class Thermo>
inline void Foam::tabulated2Transport<Thermo>::operator-=
(
    const tabulated2Transport<Thermo>& tt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(tt);
    /*
    molr1 /= this->nMoles();
    scalar molr2 = st.nMoles()/this->nMoles();

    mu_ = molr1*mu_ - molr2*st.mu_;
    rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    */
}


template<class Thermo>
inline void Foam::tabulated2Transport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulated2Transport<Thermo> Foam::operator+
(
    const tabulated2Transport<Thermo>& tt1,
    const tabulated2Transport<Thermo>& tt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(tt1) + static_cast<const Thermo&>(tt2)
    );

    scalar molr1 = tt1.nMoles()/t.nMoles();
    scalar molr2 = tt2.nMoles()/t.nMoles();

    return tabulated2Transport<Thermo>
    (
        t,
        tt1.mu,
	      tt1.kappa
    );
}


template<class Thermo>
inline Foam::tabulated2Transport<Thermo> Foam::operator-
(
    const tabulated2Transport<Thermo>& tt1,
    const tabulated2Transport<Thermo>& tt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(tt1) - static_cast<const Thermo&>(tt2)
    );

    scalar molr1 = tt1.nMoles()/t.nMoles();
    scalar molr2 = tt2.nMoles()/t.nMoles();

    return tabulated2Transport<Thermo>
    (
        t,
        tt1.mu,
	      tt1.kappa
    );
}


template<class Thermo>
inline Foam::tabulated2Transport<Thermo> Foam::operator*
(
    const scalar s,
    const tabulated2Transport<Thermo>& tt
)
{
    return tabulated2Transport<Thermo>
    (
        s*static_cast<const Thermo&>(tt),
        tt.mu_,
        tt.kappa_
    );
}


template<class Thermo>
inline Foam::tabulated2Transport<Thermo> Foam::operator==
(
    const tabulated2Transport<Thermo>& tt1,
    const tabulated2Transport<Thermo>& tt2
)
{
    return tt2 - tt1;
}


// ************************************************************************* //
