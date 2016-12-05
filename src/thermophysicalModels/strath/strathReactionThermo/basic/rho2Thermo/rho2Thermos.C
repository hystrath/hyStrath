/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "rho2Thermo.H"
#include "make2Thermo.H"

#include "specie.H"
#include "advancedSpecie.H"
#include "perfect2Gas.H"
#include "incompressiblePerfectGas.H"
#include "rhoConst.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"

#include "hConstThermo.H"
#include "janafThermo.H"

#include "decoupledEnergyModesThermo.H"
#include "sensible2Enthalpy.H"
#include "sensible2InternalEnergy.H"
#include "multiThermo.H"

#include "constantTransport.H"
#include "SutherlandEuckenTransport.H"
#include "BlottnerEuckenTransport.H"
#include "powerLawEuckenTransport.H"
#include "CEATransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"
//#include "tabulatedTransport.H" // NEW VINCENT TODO

#include "heRho2Thermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

/*make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    constantTransport,
    sensible2Enthalpy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    SutherlandEuckenTransport,
    sensible2Enthalpy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    BlottnerEuckenTransport,
    sensible2Enthalpy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    powerLawEuckenTransport,
    sensible2Enthalpy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    constantTransport,
    sensible2InternalEnergy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    SutherlandEuckenTransport,
    sensible2InternalEnergy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    BlottnerEuckenTransport,
    sensible2InternalEnergy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);

make2Thermo
(
    rho2Thermo,
    heRho2Thermo,
    pureMixture,
    powerLawEuckenTransport,
    sensible2InternalEnergy,
    decoupledEnergyModesThermo,
    perfect2Gas,
    advancedSpecie
);*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
