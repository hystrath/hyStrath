/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "reducedUnits.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const reducedUnits& rU)
{
    os  << nl << "Defined: " << nl
        << tab << "refLength = " << rU.refLength() << " m" << nl
        << tab << "refTime = " << rU.refTime() << " s" << nl
        << tab << "refMass = " << rU.refMass() << " kg" << nl
        << tab << "refCharge = " << rU.refCharge() << " C" << nl
        << nl << "Calculated: " << nl
        << tab << "refEnergy = " << rU.refEnergy() << " J" << nl
        << tab << "refTemp = " << rU.refTemp() << " K" << nl
        << tab << "refForce = " << rU.refForce() << " N" << nl
        << tab << "refVelocity = " << rU.refVelocity() << " m/s" << nl
        << tab << "refVolume = " << rU.refVolume() << " m^3" << nl
        << tab << "refPressure = " << rU.refPressure() << " N/m^2" << nl
        << tab << "refMassDensity = " << rU.refMassDensity() << " kg/m^3" << nl
        << tab << "refNumberDensity = " << rU.refNumberDensity() << " m^-3" << nl
        << tab << "refHeatFlux = " << rU.refHeatFlux() << " kg/s^-3" << nl
        << tab << "refAmpere = " << rU.refAmpere() << " C/s" << nl
        << nl << "Constants: " << nl
        << tab << "Boltzmann constant, kb = " << reducedUnits::kb 
               << " J/K, reduced: " << rU.kB() << nl
        << tab << "Elementary charge = " << reducedUnits::elementaryCharge 
               << " C, reduced: " << rU.epsilonCharge() << nl
        << tab << "Vacuum permittivity = " << reducedUnits::vacuumPermittivity 
               << " m^-3 kg^-1 C^2 s^2, reduced: " << rU.epsilonPermittivity() 
        << endl;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::reducedUnits&)"
    );

    return os;
}


// ************************************************************************* //
