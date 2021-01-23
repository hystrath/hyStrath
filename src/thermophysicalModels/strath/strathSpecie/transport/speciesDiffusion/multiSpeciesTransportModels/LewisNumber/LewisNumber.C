/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
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

#include "LewisNumber.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LewisNumber<ThermoType>::updateCoefficients()
{
    this->D_[0] = this->thermo_.kappa()*Le_ / this->thermo_.CpMix();

    for(int speciei=1; speciei < this->D_.size(); speciei++)
    {
        scalar factor = 1.0;
        
        if (this->thermo_.composition().isIon(speciei))
        {
            //- Ambipolar diffusion for charged particles
            //  In: G. V. Candler and I. Nompelis.
            //  "Computational Fluid Dynamics for Atmospheric Entry"
            //  RTO-AVT-VKI Lecture Series 2009-AVT-162, p. 22, 2009
            factor = 2.0;
        }
        
        this->D_[speciei] = factor*this->D_[0];
    }
    
    if (this->thermo_.composition().isIon(0))
    {
        this->D_[0] *= 2.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::LewisNumber<ThermoType>::LewisNumber
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    Fick<ThermoType>(thermo, turbulence),
    Le_
    (
        readScalar
        (
            IOdictionary::subDict("transportModels")
                .subDict("diffusionModelParameters").lookup("LewisNumber")
        )
    )
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
