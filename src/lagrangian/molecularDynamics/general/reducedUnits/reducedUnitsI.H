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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


inline bool Foam::reducedUnits::runReducedUnits() const
{
    return reducedUnits_;
}

inline bool Foam::reducedUnits::outputSIUnits() const
{
    return outputSI_;
}


inline Foam::scalar Foam::reducedUnits::refLength() const
{
    return refLength_;
}
// inline Foam::scalar& Foam::reducedUnits::refLength() 
// {
//     return refLength_;
// }


inline Foam::scalar Foam::reducedUnits::refTime() const
{
    return refTime_;
}
// inline Foam::scalar& Foam::reducedUnits::refTime() 
// {
//     return refTime_;
// }

inline Foam::scalar Foam::reducedUnits::refMass() const
{
    return refMass_;
}
// inline Foam::scalar& Foam::reducedUnits::refMass() 
// {
//     return refMass_;
// }

inline Foam::scalar Foam::reducedUnits::refCharge() const
{
    return refCharge_;
}
// inline Foam::scalar& Foam::reducedUnits::refCharge() 
// {
//     return refCharge_;
// }

inline Foam::scalar Foam::reducedUnits::refTemp() const
{
    return refTemp_;
}


inline Foam::scalar Foam::reducedUnits::refEnergy() const
{
    return refEnergy_;
}


inline Foam::scalar Foam::reducedUnits::refForce() const
{
    return refForce_;
}


inline Foam::scalar Foam::reducedUnits::refVelocity() const
{
    return refVelocity_;
}


inline Foam::scalar Foam::reducedUnits::refVolume() const
{
    return refVolume_;
}


inline Foam::scalar Foam::reducedUnits::refPressure() const
{
    return refPressure_;
}


inline Foam::scalar Foam::reducedUnits::refMassDensity() const
{
    return refMassDensity_;
}


inline Foam::scalar Foam::reducedUnits::refNumberDensity() const
{
    return refNumberDensity_;
}

inline Foam::scalar Foam::reducedUnits::refMolFlux() const
{
    return 1.0/refTime();
}

inline Foam::scalar Foam::reducedUnits::refMassFlux() const
{
    return refMass()/refTime();
}



inline Foam::scalar Foam::reducedUnits::refHeatFlux() const
{
    return refHeatFlux_;
}


inline Foam::scalar Foam::reducedUnits::refAmpere() const
{
    return refAmpere_;
}

//- calculated constants

inline Foam::scalar Foam::reducedUnits::kB() const
{
    return (kb*refTemp())/refEnergy();
}

inline Foam::scalar Foam::reducedUnits::epsilonPermittivity() const
{
    return vacuumPermittivity/
           (
                (refCharge()*refCharge())/
                (refForce()*refLength()*refLength())
           );
}

inline Foam::scalar Foam::reducedUnits::epsilonCharge() const
{
    return elementaryCharge/(refCharge_);
}


// ************************************************************************* //
