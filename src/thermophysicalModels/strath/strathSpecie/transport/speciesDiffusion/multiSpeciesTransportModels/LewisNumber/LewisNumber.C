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

#include "LewisNumber.H"
//#include <time.h>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LewisNumber<ThermoType>::updateCoefficients()
{
    // rho*Ds = Le*kappa_tr/Cp_tr - NEW VINCENT 16/05/2016
    /*std::clock_t start;
    double duration;*/

    //this->D_[0] = this->turbulence_.kappaEff()*Le_ / this->thermo_.Cp_t(); //TODO
    this->D_[0] = this->thermo_.kappatr()*Le_ / this->thermo_.Cp_t();

    /*start = std::clock();
    volScalarField kk = this->turbulence_.kappaEff();
    duration = (std::clock() - start) / double(CLOCKS_PER_SEC);
    Info << "timer kappa load: " << duration << endl;
    start = std::clock();
    volScalarField Cpt = this->thermo_.Cp_t();
    duration = (std::clock() - start) / double(CLOCKS_PER_SEC);
    Info << "timer Cpt load: " << duration << endl;*/

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
    const compressible::turbulenceModel& turbulence
)
:
    Fick<ThermoType>(thermo, turbulence),

    Le_(readScalar(IOdictionary::subDict("transportModels")
        .subDict("diffusionModelParameters").lookup("LewisNumber")))
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
