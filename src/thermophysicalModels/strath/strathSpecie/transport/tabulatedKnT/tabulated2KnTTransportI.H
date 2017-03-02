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
inline Foam::tabulated2KnTTransport<Thermo>::tabulated2KnTTransport
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
inline Foam::tabulated2KnTTransport<Thermo>::tabulated2KnTTransport
(
    const word& name,
    const tabulated2KnTTransport& tt
)
:
    Thermo(name, tt),
    mu_(tt.mu_),
    kappa_(tt.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2KnTTransport<Thermo> >
Foam::tabulated2KnTTransport<Thermo>::clone() const
{
    return autoPtr<tabulated2KnTTransport<Thermo> >
    (
        new tabulated2KnTTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2KnTTransport<Thermo> >
Foam::tabulated2KnTTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<tabulated2KnTTransport<Thermo> >
    (
        new tabulated2KnTTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::tabulated2KnTTransport<Thermo> >
Foam::tabulated2KnTTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<tabulated2KnTTransport<Thermo> >
    (
        new tabulated2KnTTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::tabulated2KnTTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_(T, p);
}


/*template<class Thermo>
inline Foam::scalar Foam::tabulated2KnTTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}*/
template<class Thermo>
inline Foam::scalar Foam::tabulated2KnTTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_(T, p);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulated2KnTTransport<Thermo>& Foam::tabulated2KnTTransport<Thermo>::operator=
(
    const tabulated2KnTTransport<Thermo>& tt
)
{
    Thermo::operator=(tt);

    mu_ = tt.mu_;
    kappa_ = tt.kappa_;

    return *this;
}


template<class Thermo>
inline void Foam::tabulated2KnTTransport<Thermo>::operator+=
(
    const tabulated2KnTTransport<Thermo>& tt
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
inline void Foam::tabulated2KnTTransport<Thermo>::operator-=
(
    const tabulated2KnTTransport<Thermo>& tt
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
inline void Foam::tabulated2KnTTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::tabulated2KnTTransport<Thermo> Foam::operator+
(
    const tabulated2KnTTransport<Thermo>& tt1,
    const tabulated2KnTTransport<Thermo>& tt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(tt1) + static_cast<const Thermo&>(tt2)
    );

    scalar molr1 = tt1.nMoles()/t.nMoles();
    scalar molr2 = tt2.nMoles()/t.nMoles();

    return tabulated2KnTTransport<Thermo>
    (
        t,
        tt1.mu,
	tt1.kappa
    );
}


template<class Thermo>
inline Foam::tabulated2KnTTransport<Thermo> Foam::operator-
(
    const tabulated2KnTTransport<Thermo>& tt1,
    const tabulated2KnTTransport<Thermo>& tt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(tt1) - static_cast<const Thermo&>(tt2)
    );

    scalar molr1 = tt1.nMoles()/t.nMoles();
    scalar molr2 = tt2.nMoles()/t.nMoles();

    return tabulated2KnTTransport<Thermo>
    (
        t,
        tt1.mu,
	      tt1.kappa
    );
}


template<class Thermo>
inline Foam::tabulated2KnTTransport<Thermo> Foam::operator*
(
    const scalar s,
    const tabulated2KnTTransport<Thermo>& tt
)
{
    return tabulated2KnTTransport<Thermo>
    (
        s*static_cast<const Thermo&>(tt),
        tt.mu_,
        tt.kappa_
    );
}


template<class Thermo>
inline Foam::tabulated2KnTTransport<Thermo> Foam::operator==
(
    const tabulated2KnTTransport<Thermo>& tt1,
    const tabulated2KnTTransport<Thermo>& tt2
)
{
    return tt2 - tt1;
}


// ************************************************************************* //
