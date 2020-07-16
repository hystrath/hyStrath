/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

Description

\*---------------------------------------------------------------------------*/

#include "pdAllConfigurations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor

//- Constructor for mdInitialise
pdAllConfigurations::pdAllConfigurations
(
//     const polyMesh& mesh,
    const IOdictionary& pdInitialiseDict,
    pdCloud& cloud
)
:
    pdInitialiseDict_(pdInitialiseDict),
    configurationList_(pdInitialiseDict_.lookup("configurations")),
    ids_(configurationList_.size()),
    configurations_(configurationList_.size())
{

    Info << nl << "Creating pd configurations: " << nl << endl;

    if(configurations_.size() > 0 )
    {
        forAll(configurations_, c)
        {
            const entry& configurationI = configurationList_[c];
            const dictionary& configurationIDict = configurationI.dict();

            configurations_[c] = autoPtr<pdConfiguration>
            (
                pdConfiguration::New(cloud, configurationIDict)
            );

            ids_[c] = c;
        }
    }
}

//- initial configuration

void pdAllConfigurations::setInitialConfig()
{
    forAll(configurations_, c)
    {
        configurations_[c]->setInitialConfiguration();
    }

}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
