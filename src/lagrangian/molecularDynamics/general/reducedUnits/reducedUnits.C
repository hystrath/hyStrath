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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::reducedUnits::kb = 1.3806504e-23;

const Foam::scalar Foam::reducedUnits::elementaryCharge = 1.602176487e-19;

const Foam::scalar Foam::reducedUnits::vacuumPermittivity = 8.854187817e-12;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reducedUnits::calcRefValues()
{
    if
    (
        refTime_ < VSMALL
     || refLength_ < VSMALL
     || refMass_ < VSMALL
     || refCharge_ < VSMALL
    )
    {
        FatalErrorIn("Foam::reducedUnits::calcRefValues() ")
            << "One or more reference values too small for floating point "
            << "calculation: "
            << "refTime_ = " << refTime_
            << ", refLength = " << refLength_
            << ", refMass = " << refMass_
            << ", refCharge = " << refCharge_
            << nl << abort(FatalError);
    }

    refEnergy_ = refLength_*refLength_*refMass_/(refTime_*refTime_);

    refTemp_ = refEnergy_ / kb;

    refForce_ = refEnergy_/refLength_;

    refVelocity_ = Foam::sqrt(refEnergy_/refMass_);

    refVolume_ = Foam::pow(refLength_,3.0);

    refPressure_ = refEnergy_/refVolume_;

    refMassDensity_ = refMass_/refVolume_;

    refNumberDensity_ = 1.0/refVolume_;

    refHeatFlux_ = refMass_/Foam::pow(refTime_,3.0);

    refAmpere_ = refCharge_/refTime_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reducedUnits::reducedUnits()
:
    reducedUnits_(true),
    outputSI_(false),
    refLength_(0.34e-9),
    refTime_(2.16059e-12),
    refMass_(1.660538782e-27),
    refCharge_(1.602176487e-19),
    refEnergy_(0.0),
    refTemp_(0.0),
    refForce_(0.0),
    refVelocity_(0.0),
    refVolume_(0.0),
    refPressure_(0.0),
    refMassDensity_(0.0),
    refNumberDensity_(0.0),
    refHeatFlux_(0.0),
    refAmpere_(0.0)
{

    calcRefValues();
}


Foam::reducedUnits::reducedUnits
(
    label unity
)
:
    reducedUnits_(true),
    outputSI_(false),
    refLength_(1),
    refTime_(1),
    refMass_(1),
    refCharge_(1),
    refEnergy_(0.0),
    refTemp_(0.0),
    refForce_(0.0),
    refVelocity_(0.0),
    refVolume_(0.0),
    refPressure_(0.0),
    refMassDensity_(0.0),
    refNumberDensity_(0.0),
    refHeatFlux_(0.0),
    refAmpere_(0.0)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits
(
    scalar refLength,
    scalar refTime,
    scalar refMass,
    scalar refCharge
)
:
    reducedUnits_(true),
    outputSI_(false),
    refLength_(refLength),
    refTime_(refTime),
    refMass_(refMass),
    refCharge_(refCharge),
    refEnergy_(0.0),
    refTemp_(0.0),
    refForce_(0.0),
    refVelocity_(0.0),
    refVolume_(0.0),
    refPressure_(0.0),
    refMassDensity_(0.0),
    refNumberDensity_(0.0),
    refHeatFlux_(0.0),
    refAmpere_(0.0)
{
    calcRefValues();
}


Foam::reducedUnits::reducedUnits(const IOdictionary& reducedUnitsDict)
:
    reducedUnits_(true),
    outputSI_(false),
    refLength_(0.0),
    refTime_(0.0),
    refMass_(0.0),
    refCharge_(0.0),
    refEnergy_(0.0),
    refTemp_(0.0),
    refForce_(0.0),
    refVelocity_(0.0),
    refVolume_(0.0),
    refPressure_(0.0),
    refMassDensity_(0.0),
    refNumberDensity_(0.0),
    refHeatFlux_(0.0),
    refAmpere_(0.0)
{
    setRefValues(reducedUnitsDict);
}


Foam::reducedUnits::reducedUnits
(
		Time& runTime,
		const polyMesh& mesh
)
:
    reducedUnits_(true),
    outputSI_(false),
    refLength_(0.0),
    refTime_(0.0),
    refMass_(0.0),
    refCharge_(0.0),
    refEnergy_(0.0),
    refTemp_(0.0),
    refForce_(0.0),
    refVelocity_(0.0),
    refVolume_(0.0),
    refPressure_(0.0),
    refMassDensity_(0.0),
    refNumberDensity_(0.0),
    refHeatFlux_(0.0),
    refAmpere_(0.0)
{
    IOobject reducedUnitsDictIOobject
    (
        "reducedUnitsDict",
        runTime.system(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    
    if (true/*reducedUnitsDictIOobject.headerOk()*/)
    {
        Info << nl
            << "Reading reference quantities from reducedUnitsDict file." << endl;
    
        IOdictionary reducedUnitsDict(reducedUnitsDictIOobject);

        setRefValues(reducedUnitsDict);

        Info << *this << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reducedUnits::~reducedUnits()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::reducedUnits::setRefValues
(
    scalar refLength,
    scalar refTime,
    scalar refMass,
    scalar refCharge
)
{
    refLength_ = refLength;

    refTime_ = refTime;

    refMass_ = refMass;

    refCharge_ = refCharge;

    calcRefValues();
}


void Foam::reducedUnits::setRefValues
(
    const IOdictionary& reducedUnitsDict
)
{
    reducedUnits_ = Switch(reducedUnitsDict.lookup("reducedUnits"));

    outputSI_ = Switch(reducedUnitsDict.lookup("outputSI"));

    refLength_ = readScalar(reducedUnitsDict.lookup("refLength"));

    refTime_ = readScalar(reducedUnitsDict.lookup("refTime"));

    refMass_  = readScalar(reducedUnitsDict.lookup("refMass"));

    refCharge_  = readScalar(reducedUnitsDict.lookup("refCharge"));

    calcRefValues();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::reducedUnits::operator=(const reducedUnits& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::reducedUnits::operator=(const Foam::reducedUnits&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// ************************************************************************* //
