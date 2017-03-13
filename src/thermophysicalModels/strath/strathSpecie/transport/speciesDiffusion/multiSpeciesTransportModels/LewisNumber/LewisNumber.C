/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "LewisNumber.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LewisNumber<ThermoType>::updateCoefficients()
{
    // rho*Ds = Le*kappa_tr/Cp_tr - NEW VINCENT 16/05/2016
    this->D_[0] = Le_*this->turbulence_.kappaEff() / this->thermo_.Cp_t();
    
    if(this->thermo_.composition().particleType(0) == 3)
    {
        this->D_[0] *= 2.0; // Ambipolar diffusion for charged particles
    }
    
    for(int speciei=1; speciei < this->D_.size(); speciei++)
    {
        this->D_[speciei] = this->D_[0];
        
        if(this->thermo_.composition().particleType(speciei) == 3)
        {
            this->D_[speciei] *= 2.0; // Ambipolar diffusion for charged particles
        }
    }
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::LewisNumber<ThermoType>::LewisNumber
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    Fick<ThermoType>(thermo, turbulence),
    
    Le_(readScalar(IOdictionary::subDict("transportModels").lookup("LewisNumber"))) 
{}
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
