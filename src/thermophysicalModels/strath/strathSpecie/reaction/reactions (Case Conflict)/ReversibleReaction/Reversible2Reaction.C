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

#include "Reversible2Reaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::Reversible2Reaction<ReactionType, ReactionThermo, ReactionRate>::
Reversible2Reaction
(
    const ReactionType<ReactionThermo>& reaction,
    const ReactionRate& k,
    const DynamicList<scalar> ni, // NEW VINCENT 09/02/2016
    const DynamicList<scalar> A0,
    const DynamicList<scalar> A1,
    const DynamicList<scalar> A2,
    const DynamicList<scalar> A3,
    const DynamicList<scalar> A4
)
:
    ReactionType<ReactionThermo>(reaction),
    k_(k),
    ni_(ni), // NEW VINCENT 09/02/2016
    A0_(A0),
    A1_(A1),
    A2_(A2),
    A3_(A3),
    A4_(A4)
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::Reversible2Reaction<ReactionType, ReactionThermo, ReactionRate>::
Reversible2Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, is),
    k_(species, is)
{

    forAll(ni_, i)
    {
        is >> ni_[i];
    }
    
    forAll(A0_, i)
    {
        is >> A0_[i];
    }
    
    forAll(A1_, i)
    {
        is >> A1_[i];
    }
    
    forAll(A2_, i)
    {
        is >> A2_[i];
    }
    
    forAll(A3_, i)
    {
        is >> A3_[i];
    }
    
    forAll(A4_, i)
    {
        is >> A4_[i];
    }

}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::Reversible2Reaction<ReactionType, ReactionThermo, ReactionRate>::
Reversible2Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict),
    ni_(dict.lookup("ni")), // NEW VINCENT 09/02/2016
    A0_(dict.lookup("A0")),
    A1_(dict.lookup("A1")),
    A2_(dict.lookup("A2")),
    A3_(dict.lookup("A3")),
    A4_(dict.lookup("A4"))
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::Reversible2Reaction<ReactionType, ReactionThermo, ReactionRate>::
Reversible2Reaction
(
    const Reversible2Reaction<ReactionType, ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    ReactionType<ReactionThermo>(rr, species),
    k_(rr.k_),
    ni_(rr.ni_), // NEW VINCENT 09/02/2016
    A0_(rr.A0_),
    A1_(rr.A1_),
    A2_(rr.A2_),
    A3_(rr.A3_),
    A4_(rr.A4_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::Reversible2Reaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return k_(p, T, c);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::Reversible2Reaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    //return kfwd/this->Kc(p, T); // ORIGINAL FORMULATION, DELETED VINCENT 09/02/2016
    return kfwd/Keq(p, T); // NEW VINCENT 09/02/2016
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::Reversible2Reaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return kr(kf(p, T, c), p, T, c);
}


// NEW VINCENT 09/02/2016 *****************************************************
template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::Reversible2Reaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::Keq
(
    const scalar p,
    const scalar T
) const
{
    scalar nDMix = p/(constant::physicoChemical::k.value()*T);
    const label nLow(findBounds(nDMix)), nHigh(nLow+1);
    
    const scalar A0Mix = boundedLinearInterpolation(nDMix, ni_[nLow], ni_[nHigh], A0_[nLow], A0_[nHigh]);
    const scalar A1Mix = boundedLinearInterpolation(nDMix, ni_[nLow], ni_[nHigh], A1_[nLow], A1_[nHigh]);
    const scalar A2Mix = boundedLinearInterpolation(nDMix, ni_[nLow], ni_[nHigh], A2_[nLow], A2_[nHigh]);
    const scalar A3Mix = boundedLinearInterpolation(nDMix, ni_[nLow], ni_[nHigh], A3_[nLow], A3_[nHigh]);
    const scalar A4Mix = boundedLinearInterpolation(nDMix, ni_[nLow], ni_[nHigh], A4_[nLow], A4_[nHigh]);
    
    const scalar Z = 1.0e4/T;
    
    //Info << "Keq: " << exp(A0Mix*T/1.0e4 + A1Mix + A2Mix*log(1.0e4/T) + A3Mix*1.0e4/T + A4Mix*sqr(1.0e4/T)) << endl;
    return exp(A0Mix/Z + A1Mix + A2Mix*log(Z) + A3Mix*Z + A4Mix*sqr(Z));
}
// END NEW VINCENT 09/02/2016 *************************************************


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::Reversible2Reaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::write
(
    Ostream& os
) const
{
    Reaction2<ReactionThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
