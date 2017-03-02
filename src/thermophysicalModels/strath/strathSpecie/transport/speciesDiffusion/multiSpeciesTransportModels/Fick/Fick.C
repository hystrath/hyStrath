/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "Fick.H"
#include "fvm.H"
#include <ctime>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::Fick<ThermoType>::updateCoefficients()
{     
    DijModel_().update();

    forAll(species(), speciei)
    {
        volScalarField tmpSum = 0 / Dij(0,0);

        forAll(species(), speciej)
        {
            if (speciej != speciei and thermo_.composition().particleType(speciej) != 0)
            {     
                tmpSum += thermo_.composition().X(speciej) / Dij(speciei, speciej);
            }
        }
        
        volScalarField& Xi = thermo_.composition().X(speciei);
        
        D_[speciei] = thermo_.rho()*(1.0 - Xi) 
            / (tmpSum + dimensionedScalar("VSMALL", dimTime/dimArea, Foam::VSMALL));   

        
        forAll(D_[speciei], celli)
        {
            if (1.0 - Xi[celli] < miniXs_)
            {
                D_[speciei][celli] = 0;
            }
        }
        
        forAll(D_[speciei].boundaryField(), patchi)
        {
            forAll(D_[speciei].boundaryField()[patchi], facei)
            {
                if (1.0 - Xi.boundaryField()[patchi][facei] < miniXs_)
                {
                    D_[speciei].boundaryField()[patchi][facei] = 0;
                }  
            }
        }
    }
} 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Fick<ThermoType>::Fick
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    multiSpeciesTransportModel(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),
    
    miniXs_(1.0e-4)
{    
    D_.setSize(species().size());
    
    forAll(species(), speciei)
    {
        D_.set
        (
            speciei, 
            new volScalarField
            (
                IOobject
                (
                    "rhoD_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("D", dimMass/dimLength/dimTime, 0.0)
            )
        );
    } 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::Fick<ThermoType>::correct()
{
    updateCoefficients();

    forAll(species(), speciei)
    {
        calculateJ(speciei);
    }
    
    calculateSumDiffusiveFluxes();
}

    
template<class ThermoType>
bool Foam::Fick<ThermoType>::read()
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
