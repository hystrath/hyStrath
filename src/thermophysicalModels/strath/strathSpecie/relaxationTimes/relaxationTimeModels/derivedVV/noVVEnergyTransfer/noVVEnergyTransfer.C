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

#include "noVVEnergyTransfer.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::noVVEnergyTransfer<ThermoType>::updateCoefficients()
{     
    forAll(species(), i)
    {
        tauVV_[i] = dimensionedScalar("GREAT", dimTime, Foam::GREAT);
    }
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::noVVEnergyTransfer<ThermoType>::noVVEnergyTransfer
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    relaxationTimeModelVV(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{    
    tauVV_.setSize(solvedVibEqSpecies().size()); // NEW VINCENT 06/08/2016
    
    forAll(tauVV_, speciei)
    {
        tauVV_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "tauVV_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("tauVV", dimTime, 0.0)
            )
        );
    } 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::noVVEnergyTransfer<ThermoType>::correct()
{}  
    

template<class ThermoType>
bool Foam::noVVEnergyTransfer<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}
   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
