/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "decoupledEnergyModesThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::decoupledEnergyModesThermo<EquationOfState>::checkInputData() const
{
    if (Tlow_ >= Thigh_)
    {
        FatalErrorIn("decoupledEnergyModesThermo<EquationOfState>::check()")
            << "Tlow(" << Tlow_ << ") >= Thigh(" << Thigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::decoupledEnergyModesThermo<EquationOfState>::decoupledEnergyModesThermo(Istream& is)
:
    EquationOfState(is),
    Tlow_(readScalar(is)),
    Thigh_(readScalar(is))
{
    checkInputData();

    forAll(decoupledCvCoeffs_, i)
    {
        is >> decoupledCvCoeffs_[i];
    }
    
    forAll(vibrationalList_, i)
    {
        is >> vibrationalList_[i];
    }
    
    forAll(electronicList_, i)
    {
        is >> electronicList_[i];
    }

    // Check state of Istream
    is.check("decoupledEnergyModesThermo::decoupledEnergyModesThermo(Istream& is)");
}


template<class EquationOfState>
Foam::decoupledEnergyModesThermo<EquationOfState>::decoupledEnergyModesThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tlow_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("Tlow", Foam::SMALL)),
    Thigh_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("Thigh", Foam::GREAT)),
    decoupledCvCoeffs_(dict.subDict("thermodynamics").lookup("decoupledCvCoeffs")),
    vibrationalList_(dict.subDict("thermodynamics").lookup("vibrationalList")),
    electronicList_(dict.subDict("thermodynamics").lookup("electronicList"))
{
    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::decoupledEnergyModesThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Tlow", Tlow_);
    dict.add("Thigh", Thigh_);
    dict.add("decoupledCvCoeffs", decoupledCvCoeffs_);
    dict.add("vibrationalList", vibrationalList_);
    dict.add("electronicList", electronicList_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const decoupledEnergyModesThermo<EquationOfState>& dem
)
{
    os  << static_cast<const EquationOfState&>(dem) << nl
        << "    " << dem.Tlow_
        << tab << dem.Thigh_;

    os << nl << "    ";

    forAll(dem.decoupledCvCoeffs_, i)
    {
        os << dem.decoupledCvCoeffs_[i] << ' ';
    }
    
    os << nl << "    ";

    forAll(dem.vibrationalList_, i)
    {
        os << dem.vibrationalList_[i] << ' ';
    }
    
    os << nl << "    ";

    forAll(dem.electronicList_, i)
    {
        os << dem.electronicList_[i] << ' ';
    }

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const decoupledEnergyModesThermo<EquationOfState>& dem)"
    );

    return os;
}


// ************************************************************************* //
