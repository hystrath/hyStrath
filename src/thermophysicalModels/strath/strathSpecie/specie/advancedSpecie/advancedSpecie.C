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

#include "advancedSpecie.H"
#include "IOstreams.H"
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

//- Universal gas constant (default in [J/(kmol K)])
const Foam::scalar Foam::advancedSpecie::RR = constant::physicoChemical::R.value()*1000;

//- Standard pressure (default in [Pa])
const Foam::scalar Foam::advancedSpecie::Pstd = constant::standard::Pstd.value();

//- Standard temperature (default in [K])
const Foam::scalar Foam::advancedSpecie::Tstd = constant::standard::Tstd.value();


namespace Foam
{
    defineTypeNameAndDebug(advancedSpecie, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::advancedSpecie::advancedSpecie(Istream& is)
:
    name_(is),
    nMoles_(readScalar(is)),
    molWeight_(readScalar(is)),
    particleType_(readScalar(is)),
    particleCharge_(readScalar(is)),
    diameter_(readScalar(is)),
    omega_(readScalar(is)),
    dissociationPotential_(readScalar(is)),
    noVibrationalTemp_(readScalar(is)),
    noElectronicLevels_(readScalar(is)),
    iHat_(readScalar(is))
{
    is.check("advancedSpecie::advancedSpecie(Istream& is)");
    
    forAll(vibrationalList_, i)
    {
        is >> vibrationalList_[i];
    }
}


Foam::advancedSpecie::advancedSpecie(const dictionary& dict)
:
    name_(dict.dictName()),
    nMoles_(readScalar(dict.subDict("specie").lookup("nMoles"))),
    molWeight_(readScalar(dict.subDict("specie").lookup("molWeight"))),
    particleType_(readScalar(dict.subDict("specie").lookup("particleType"))),
    particleCharge_(dict.subDict("specie").lookupOrDefault<scalar>("charge", 0)),
    diameter_(readScalar(dict.subDict("specie").lookup("diameter"))),
    omega_(readScalar(dict.subDict("specie").lookup("omega"))),
    vibrationalList_(dict.subDict("thermodynamics").lookup("vibrationalList")),
    dissociationPotential_(readScalar(dict.subDict("specie").lookup("dissocEnergy"))),
    noVibrationalTemp_(readScalar(dict.subDict("specie").lookup("noVibTemp"))),
    noElectronicLevels_(readScalar(dict.subDict("specie").lookup("noElecLevels"))),
    iHat_(dict.subDict("specie").lookupOrDefault<scalar>("iHat", 0))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::advancedSpecie::write(Ostream& os) const
{
    dictionary dict("specie");
    dictionary dict2("thermodynamics");
    dict.add("nMoles", nMoles_);
    dict.add("molWeight", molWeight_);
    dict.add("particleType", particleType_);
    dict.add("charge", particleCharge_);
    dict.add("diameter", diameter_);
    dict.add("omega", omega_);
    dict2.add("vibrationalList", vibrationalList_);
    dict.add("dissociationPotential", dissociationPotential_);
    dict.add("noVibrationalTemp", noVibrationalTemp_);
    dict.add("noElectronicLevels", noElectronicLevels_);
    dict.add("iHat", iHat_);
    os  << indent << dict.dictName() << dict;
    os  << indent << dict2.dictName() << dict2;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const advancedSpecie& as)
{
    os  << as.name_ << tab
        << as.nMoles_ << tab
        << as.molWeight_ << tab
        << as.particleType_ << tab
        << as.particleCharge_ << tab
        << as.diameter_ << tab
        << as.omega_ << tab;
        
        forAll(as.vibrationalList_, i)
        {
            os << as.vibrationalList_[i] << tab;
        }
        
    os  << as.dissociationPotential_ << tab
        << as.noVibrationalTemp_ << tab
        << as.noElectronicLevels_ << tab
        << as.iHat_ << tab;

    os.check("Ostream& operator<<(Ostream& os, const advancedSpecie& as)");
    return os;
}


// ************************************************************************* //
