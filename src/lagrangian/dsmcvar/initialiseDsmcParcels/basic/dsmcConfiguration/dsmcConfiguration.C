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

#include "dsmcConfiguration.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcConfiguration, 0);

defineRunTimeSelectionTable(dsmcConfiguration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcConfiguration::dsmcConfiguration
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    dsmcInitialiseDict_(dict),
    rndGen_(cloud.rndGen()),
    nParcelsAdded_(0)
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcConfiguration> dsmcConfiguration::New
(

    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcConfigurationName
    (
        dict.lookup("type")
    );

    Info<< "Selecting dsmcConfiguration "
         << dsmcConfigurationName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcConfigurationName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcConfiguration::New(const dictionary&) : " << endl
            << "    unknown dsmcConfiguration type "
            << dsmcConfigurationName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcConfiguration>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcConfiguration::~dsmcConfiguration()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// vector dsmcConfiguration::equipartitionLinearVelocity
// (
//     scalar temperature,
//     scalar mass
// )
// {
//     return sqrt(molCloud_.redUnits().kB()*temperature/mass)*vector
//     (
//         rndGen_.GaussNormal(),
//         rndGen_.GaussNormal(),
//         rndGen_.GaussNormal()
//     );
// }

// const word& dsmcConfiguration::name() const
// {
//     return name_;
// }

const label& dsmcConfiguration::nParcelsAdded() const
{
    return nParcelsAdded_;
}


} // End namespace Foam

// ************************************************************************* //
