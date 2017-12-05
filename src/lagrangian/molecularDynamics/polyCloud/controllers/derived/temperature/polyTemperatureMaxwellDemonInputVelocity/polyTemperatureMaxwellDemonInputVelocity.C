/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "polyTemperatureMaxwellDemonInputVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyTemperatureMaxwellDemonInputVelocity, 0);

addToRunTimeSelectionTable(polyStateController, polyTemperatureMaxwellDemonInputVelocity, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyTemperatureMaxwellDemonInputVelocity::polyTemperatureMaxwellDemonInputVelocity
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    tauT_(readScalar(propsDict_.lookup("tauT"))),
    p_(0.0),
    molIds_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;

    temperature_ = readScalar(propsDict_.lookup("temperature"));

    vector direction = propsDict_.lookup("velocityDirection");
    direction /= mag(direction);
    
    const word distributionName = propsDict_.lookup("velocityDistributionName");

    fileName timePath(time_.time().system()/distributionName);

    IFstream file(timePath);

    List< Pair<scalar> > velocity;

    if (file.good())
    {
        file >> velocity;
    }
    else
    {
        FatalErrorIn
        (
            "polyTemperatureMaxwellDemonInputVelocity"
        )
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
    
    label nBins = velocity.size();
    nBins_ = nBins;
    
    velocities_.setSize(nBins, vector::zero);
    
    forAll(velocity, i)
    {
        velocities_[i]=velocity[i].second()*direction;
    }
    
    binWidth_ = velocity[1].first()-velocity[0].first();
    
    Info << "Binwidth = " << binWidth_ << endl;
    
    zeroPoint_ = readScalar(propsDict_.lookup("zeroPoint"));
    
    p_ = 1.0 - exp(-(time_.deltaT().value()/tauT_));

    Info << "probability of collisions:  " << p_ << endl;

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    Info << "mol ids = " << molIds_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyTemperatureMaxwellDemonInputVelocity::~polyTemperatureMaxwellDemonInputVelocity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyTemperatureMaxwellDemonInputVelocity::initialConfiguration()
{}

void polyTemperatureMaxwellDemonInputVelocity::controlBeforeVelocityI()
{}

void polyTemperatureMaxwellDemonInputVelocity::controlBeforeMove()
{}

void polyTemperatureMaxwellDemonInputVelocity::controlBeforeForces()
{}

void polyTemperatureMaxwellDemonInputVelocity::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyTemperatureMaxwellDemonInputVelocity::controlAfterForces()
{}

void polyTemperatureMaxwellDemonInputVelocity::controlAfterVelocityII()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            vector velocity = vector::zero;
            
            scalar y = abs(zeroPoint_ - mol().position().y());
            
            label index = label(y/binWidth_);
            
//             Info << "index = " << index << endl;
            
            if( (index < nBins_) && (index >= 0) )
            {
                velocity = velocities_[index];
            }
            else
            {
                Info << "WARNING out of bounds, y = " << y << endl;
            }
            
            
            if(molCloud_.rndGen().sample01<scalar>() <= p_)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());

                scalar sigma = sqrt (temperature_*molCloud_.redUnits().kB() / massI);

                vector molVel
                (
                    sigma*molCloud_.rndGen().GaussNormalMD<scalar>(),
                    sigma*molCloud_.rndGen().GaussNormalMD<scalar>(),
                    sigma*molCloud_.rndGen().GaussNormalMD<scalar>()
                );

                mol().v() = velocity + molVel;
            }
        }
    }  
}


void polyTemperatureMaxwellDemonInputVelocity::calculateProperties()
{}

void polyTemperatureMaxwellDemonInputVelocity::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}

void polyTemperatureMaxwellDemonInputVelocity::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
/*
    if (propsDict_.found("tauT"))
    {
        tauT_ = readScalar(propsDict_.lookup("tauT"));
    }

    p_ = 1.0 - exp(-(time_.deltaT().value()/tauT_));

    Info << "probability of collisions:  " << p_ << endl;

    if (readStateFromFile_)
    {
        velocity_ = propsDict_.lookup("velocity");

        forAll(velocities_, c)
        {
            velocities_[c] = velocity_;
        }

        temperature_ = readScalar(propsDict_.lookup("temperature"));

        forAll(temperatures_, c)
        {
            temperatures_[c] = temperature_;
        }
    }*/
}

} // End namespace Foam

// ************************************************************************* //
