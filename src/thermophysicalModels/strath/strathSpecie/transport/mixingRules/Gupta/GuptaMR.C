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

#include "GuptaMR.H"
#include "fvm.H"

#include<sstream>

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

// Avogadro's number
template<class ThermoType>
const Foam::scalar Foam::GuptaMR<ThermoType>::NA = Foam::constant::physicoChemical::NA.value();

// Boltzmann's constant
template<class ThermoType>
const Foam::scalar Foam::GuptaMR<ThermoType>::kB = Foam::constant::physicoChemical::k.value();

//- Universal gas constant (in [J/(mol K)])
template<class ThermoType>
const Foam::scalar Foam::GuptaMR<ThermoType>::Ru = Foam::constant::physicoChemical::R.value();

//- Mathematical constant Pi
template<class ThermoType>
const Foam::scalar Foam::GuptaMR<ThermoType>::pi = Foam::constant::mathematical::pi;

//- Fundamental electric charge in CGS units
template<class ThermoType>
const Foam::scalar Foam::GuptaMR<ThermoType>::eCGS = 4.8032e-10;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template <class ThermoType>
Foam::string Foam::GuptaMR<ThermoType>::numberToString(int number)
{
   std::ostringstream ss;
   ss << number;
   return ss.str();
}
  

template<class ThermoType>
void Foam::GuptaMR<ThermoType>::piOmegaNeutralInit(label no, label i, label j)
{
    // NOT AVAILABLE IN BETA RELEASE
}


template<class ThermoType>
void Foam::GuptaMR<ThermoType>::piOmegaNonNeutralsInit(const word& attractionType, label no)
{
    // NOT AVAILABLE IN BETA RELEASE
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::GuptaMR<ThermoType>::GuptaMR
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel2& turbulence
)
:
    mixingRule(thermo, turbulence),
    
    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    ),
    
    dict_
    (
        thermo.transportDictionary()
    )
{
    
    label noUnchargedParticles = 0;
    forAll(species(), speciei)
    {
        if(speciesThermo_[speciei].particleCharge() == 0)
        {
            noUnchargedParticles += 1;
        }
    }
    
    forAll(piOmegaNeutral_, k)
    {
        piOmegaNeutral_[k].setSize(noUnchargedParticles);
        
        forAll(piOmegaNeutral_[k], speciei)
        { 
            piOmegaNeutral_[k].set
            (
                speciei, 
                new PtrList<FixedList<scalar, 4> >(noUnchargedParticles)
            );
        }
    }
    
    forAll(piOmegaNeutral_, k)
    {
        forAll(piOmegaNeutral_[k], speciei)
        { 
            forAll(piOmegaNeutral_[k][speciei], speciej)
            {
                piOmegaNeutral_[k][speciei].set
                (
                    speciej, 
                    new FixedList<scalar, 4>
                );
            }
        }
    }
    
    for(int speciei = 0; speciei < noUnchargedParticles; speciei++)
    { 
        for(int speciej = 0; speciej < noUnchargedParticles; speciej++)
        {
            for(int k = 0; k <= 1; k++)
            {
                piOmegaNeutralInit(k, speciei, speciej);
            }
        }
    }
    
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::GuptaMR<ThermoType>::correct()
{
    // NOT AVAILABLE IN BETA RELEASE
} 


template<class ThermoType>
void Foam::GuptaMR<ThermoType>::write()
{
    // NOT AVAILABLE IN BETA RELEASE    
} 


template<class ThermoType>
bool Foam::GuptaMR<ThermoType>::read()
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
