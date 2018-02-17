/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

#include "directionalPressureGradientExplicitSource.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::directionalPressureGradientExplicitSource::writeData
(
    Ostream& os
) const
{
    NotImplemented;
}


bool Foam::fv::directionalPressureGradientExplicitSource::read
(
    const dictionary& dict
)
{
    const dictionary coeffs(dict.subDict(typeName + "Coeffs"));

    relaxationFactor_ =
        coeffs.lookupOrDefault<scalar>("relaxationFactor", 0.3);

    coeffs.lookup("flowDir") >> flowDir_;
    flowDir_ /= mag(flowDir_);

    if (model_ == pConstant)
    {
        coeffs.lookup("pressureDrop") >> pressureDrop_;
    }
    else if (model_ == pDarcyForchheimer)
    {
        coeffs.lookup("D") >> D_;
        coeffs.lookup("I") >> I_;
        coeffs.lookup("length") >> length_;
    }

    return false;
}


// ************************************************************************* //
