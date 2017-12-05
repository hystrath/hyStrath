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

#include "polyConfiguration.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyConfiguration, 0);
defineRunTimeSelectionTable(polyConfiguration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyConfiguration::polyConfiguration
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
    mdInitialiseDict_(dict),
    nMolsAdded_(0)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyConfiguration> polyConfiguration::New
(

    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyConfigurationName
    (
        dict.lookup("type")
    );

    Info<< "Selecting polyConfiguration "
         << polyConfigurationName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyConfigurationName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyConfiguration::New(const dictionary&) : " << endl
            << "    unknown polyConfiguration type "
            << polyConfigurationName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyConfiguration>
	(
		cstrIter()(molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyConfiguration::~polyConfiguration()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector polyConfiguration::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return sqrt(molCloud_.redUnits().kB()*temperature/mass)*vector
    (
        molCloud_.rndGen().GaussNormalMD<scalar>(),
        molCloud_.rndGen().GaussNormalMD<scalar>(),
        molCloud_.rndGen().GaussNormalMD<scalar>()
    );
}

vector polyConfiguration::equipartitionAngularMomentum
(
    scalar temperature,
    label id
)
{
    scalar sqrtKbT = sqrt(molCloud_.redUnits().kB()*temperature);

    if (molCloud_.cP().linearMolecule(id))
    {
        return sqrtKbT*vector
        (
            0.0,
            sqrt(molCloud_.cP().momentOfInertia(id).yy())*molCloud_.rndGen().GaussNormalMD<scalar>(),
            sqrt(molCloud_.cP().momentOfInertia(id).zz())*molCloud_.rndGen().GaussNormalMD<scalar>()
        );
    }
    else
    {
        return sqrtKbT*vector
        (
            sqrt(molCloud_.cP().momentOfInertia(id).xx())*molCloud_.rndGen().GaussNormalMD<scalar>(),
            sqrt(molCloud_.cP().momentOfInertia(id).yy())*molCloud_.rndGen().GaussNormalMD<scalar>(),
            sqrt(molCloud_.cP().momentOfInertia(id).zz())*molCloud_.rndGen().GaussNormalMD<scalar>()
        );
    }
}


void polyConfiguration::insertMolecule
(
    const point& position,
    const label cell,
    const label tetFace,
    const label tetPt, 
    const label& id, 
    const bool& tethered,
    const bool& frozen,
    const scalar& temperature,
    const vector& bulkVelocity
)
{
    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_TETHERED;
    }

    if (frozen)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_FROZEN;
    }

//     const polyMolecule::constantProperties& cP = molCloud_.constProps(id);

    vector v = equipartitionLinearVelocity(temperature, molCloud_.cP().mass(id));

    v += bulkVelocity;

    vector pi = vector::zero;

    tensor Q = I;

    if (!molCloud_.cP().pointMolecule(id))
    {
        pi = equipartitionAngularMomentum(temperature, id);
        scalar phi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
        scalar theta(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
        scalar psi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);

        Q = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );
    }

    molCloud_.createMolecule
    (
        position,
        cell,
        tetFace,
        tetPt,     
        Q,
        v,
        vector::zero,
        pi,
        vector::zero,
        specialPosition,
        special,
        id,
        1.0,
        molCloud_.getTrackingNumber()
    );
}

void polyConfiguration::insertMolecule
(
    const point& position,
    const label& id,
    const bool& tethered,
    const bool& frozen,
    const scalar& temperature,
    const vector& bulkVelocity
)
{
    label cell = -1;
    label tetFace = -1;
    label tetPt = -1;

    mesh_.findCellFacePt
    (
        position,
        cell,
        tetFace,
        tetPt
    );    
    
    if(cell != -1)
    {
        point specialPosition(vector::zero);

        label special = 0;

        if (tethered)
        {
            specialPosition = position;

            special = polyMolecule::SPECIAL_TETHERED;
        }

        if (frozen)
        {
            specialPosition = position;

            special = polyMolecule::SPECIAL_FROZEN;
        }

//         const polyMolecule::constantProperties& cP = molCloud_.constProps(id);

        vector v = equipartitionLinearVelocity(temperature, molCloud_.cP().mass(id));

        v += bulkVelocity;

        vector pi = vector::zero;

        tensor Q = I;

        if (!molCloud_.cP().pointMolecule(id))
        {
            pi = equipartitionAngularMomentum(temperature, id);
            scalar phi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
            scalar theta(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);
            scalar psi(molCloud_.rndGen().sample01<scalar>()*constant::mathematical::twoPi);

            Q = tensor
            (
                cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                sin(psi)*sin(theta),
                - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                cos(psi)*sin(theta),
                sin(theta)*sin(phi),
                - sin(theta)*cos(phi),
                cos(theta)
            );
        }

        molCloud_.createMolecule
        (
            position,
            cell,
            tetFace,
            tetPt,     
            Q,
            v,
            vector::zero,
            pi,
            vector::zero,
            specialPosition,
            special,
            id,
            1.0,
            molCloud_.getTrackingNumber()
        );
    }
    else
    {
        Info << "WARNING. Molecule not inserted since position out of mesh = "
            << position
            << endl;
    }
}

void polyConfiguration::deleteMolecule
(
    polyMolecule& mol
)
{
    //- remove polyMolecule from cloud
    molCloud_.deleteParticle(mol);
}


const label& polyConfiguration::mols() const
{
    return nMolsAdded_;
}

} // End namespace Foam

// ************************************************************************* //
