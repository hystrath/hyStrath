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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ImpactIonisationArrhenius2ReactionRate::ImpactIonisationArrhenius2ReactionRate
(
    const DynamicList<scalar> Tbound,
    const DynamicList<scalar> A,
    const DynamicList<scalar> beta,
    const DynamicList<scalar> Ta,
    const DynamicList<scalar> C,
    const scalar smoothingTemperatureRange
)
:
    Tbound_(Tbound),
    A_(A),
    beta_(beta),
    Ta_(Ta),
    C_(C),
    smoothingTemperatureRange_(smoothingTemperatureRange)
{}


inline Foam::ImpactIonisationArrhenius2ReactionRate::ImpactIonisationArrhenius2ReactionRate
(
    const speciesTable&,
    Istream& is
)
:
    smoothingTemperatureRange_(readScalar(is))
{
    forAll(Tbound_, i)
    {
        is >> Tbound_[i];
    }
    
    forAll(A_, i)
    {
        is >> A_[i];
    }
    
    forAll(beta_, i)
    {
        is >> beta_[i];
    }
    
    forAll(Ta_, i)
    {
        is >> Ta_[i];
    }
    
    forAll(C_, i)
    {
        is >> C_[i];
    }
    
    is.readEnd("ImpactIonisationArrhenius2ReactionRate(Istream&)");
}


inline Foam::ImpactIonisationArrhenius2ReactionRate::ImpactIonisationArrhenius2ReactionRate
(
    const speciesTable&,
    const dictionary& dict
)
:
    Tbound_(dict.lookup("Tbound")),
    A_(dict.lookup("A")),
    beta_(dict.lookup("beta")),
    Ta_(dict.lookup("Ta")),
    C_(dict.lookup("C")),
    smoothingTemperatureRange_(readScalar(dict.lookup("smoothingTemperatureRange")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::ImpactIonisationArrhenius2ReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField&
) const
{
    // NEW VINCENT 03/06/2016 *************************************************
    const label temperatureRange(findBounds(T));
    
    scalar A_range = A_[temperatureRange];
    scalar beta_range = beta_[temperatureRange];
    scalar Ta_range = Ta_[temperatureRange];
    scalar C_range = C_[temperatureRange];
    
    scalar ak = A_range;

    if (mag(beta_range) > VSMALL)
    {
        ak *= pow(T, beta_range);
    }

    if (mag(Ta_range) > VSMALL)
    {
        ak *= exp(-Ta_range/T);
    }
    
    ak *= Ta_range/T + C_range;

    if(Tbound_.size() != 0)
    {
        if(T-Tbound_[temperatureRange-1] < smoothingTemperatureRange_ and T-Tbound_[temperatureRange-1] > 0)
        {
            A_range = A_[temperatureRange-1];
            beta_range = beta_[temperatureRange-1];
            Ta_range = Ta_[temperatureRange-1];
            C_range = C_[temperatureRange-1];
            
            scalar ak2 = A_range;

            if (mag(beta_range) > VSMALL)
            {
                ak2 *= pow(T, beta_range);
            }

            if (mag(Ta_range) > VSMALL)
            {
                ak2 *= exp(-Ta_range/T);
            }
            
            ak2 *= Ta_range/T + C_range;
            
            const scalar alpha = (Tbound_[temperatureRange-1]+smoothingTemperatureRange_ - T)/(2.0*smoothingTemperatureRange_);
            return ak*(1.0-alpha) + ak2*alpha;
        }
        else if(T-Tbound_[temperatureRange] > -smoothingTemperatureRange_ and T-Tbound_[temperatureRange] < 0)
        {
            A_range = A_[temperatureRange+1];
            beta_range = beta_[temperatureRange+1];
            Ta_range = Ta_[temperatureRange+1];
            C_range = C_[temperatureRange+1];
            
            scalar ak2 = A_range;

            if (mag(beta_range) > VSMALL)
            {
                ak2 *= pow(T, beta_range);
            }

            if (mag(Ta_range) > VSMALL)
            {
                ak2 *= exp(-Ta_range/T);
            }
            
            ak2 *= Ta_range/T + C_range;
            
            const scalar alpha = (Tbound_[temperatureRange]+smoothingTemperatureRange_ - T)/(2.0*smoothingTemperatureRange_);
            return ak*alpha + ak2*(1.0-alpha);
        }
        else
        {
            return ak;
        }
    }
    else
    {
        return ak;
    }
    // END NEW VINCENT 03/06/2016 *********************************************
}


inline void Foam::ImpactIonisationArrhenius2ReactionRate::write(Ostream& os) const
{
    os.writeKeyword("Tbound") << Tbound_ << token::END_STATEMENT << nl;
    os.writeKeyword("A") << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ta") << Ta_ << token::END_STATEMENT << nl;
    os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
    os.writeKeyword("smoothingTemperatureRange") << smoothingTemperatureRange_ << token::END_STATEMENT << nl;
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ImpactIonisationArrhenius2ReactionRate& arr
)
{
    os  << token::BEGIN_LIST
        << arr.A_ << token::SPACE << arr.beta_ << token::SPACE << arr.Ta_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
