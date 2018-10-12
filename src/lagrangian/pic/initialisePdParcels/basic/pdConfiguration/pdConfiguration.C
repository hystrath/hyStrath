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

#include "pdConfiguration.H"
#include "IFstream.H"
#include "graph.H"
#include "pdCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdConfiguration, 0);

defineRunTimeSelectionTable(pdConfiguration, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdConfiguration::pdConfiguration
(
    pdCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    pdInitialiseDict_(dict),
    rndGen_(cloud.rndGen()),
    nParcelsAdded_(0)
{

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pdConfiguration> pdConfiguration::New
(

    pdCloud& cloud,
    const dictionary& dict
)
{
    word pdConfigurationName
    (
        dict.lookup("type")
    );

    Info<< "Selecting pdConfiguration "
         << pdConfigurationName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdConfigurationName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdConfiguration::New(const dictionary&) : " << endl
            << "    unknown pdConfiguration type "
            << pdConfigurationName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pdConfiguration>
	(
		cstrIter()(cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdConfiguration::~pdConfiguration()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// vector pdConfiguration::equipartitionLinearVelocity
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

// const word& pdConfiguration::name() const
// {
//     return name_;
// }

const label& pdConfiguration::nParcelsAdded() const
{
    return nParcelsAdded_;
}


} // End namespace Foam

// ************************************************************************* //
