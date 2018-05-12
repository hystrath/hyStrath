/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "cloudInfo.H"
#include "kinematicCloud.H"
#include "dictionary.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cloudInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeTabbed(os, "nParcels");
    writeTabbed(os, "mass");
    writeTabbed(os, "Dmax");
    writeTabbed(os, "D10");
    writeTabbed(os, "D32");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::cloudInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::~cloudInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cloudInfo::read(const dictionary& dict)
{
    if (regionFunctionObject::read(dict) && logFiles::read(dict))
    {
        logFiles::resetNames(dict.lookup("clouds"));

        Info<< type() << " " << name() << ": ";
        if (writeToFile() && names().size())
        {
            Info<< "applying to clouds:" << nl;
            forAll(names(), i)
            {
                Info<< "    " << names()[i] << nl;
                writeFileHeader(files(i));
            }
            Info<< endl;
        }
        else
        {
            Info<< "no clouds to be processed" << nl << endl;
        }

        return true;
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::execute()
{
    return true;
}


bool Foam::functionObjects::cloudInfo::write()
{
    forAll(names(), i)
    {
        const word& cloudName = names()[i];

        const kinematicCloud& cloud =
            obr_.lookupObject<kinematicCloud>(cloudName);

        label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());
        scalar massInSystem =
            returnReduce(cloud.massInSystem(), sumOp<scalar>());

        scalar Dmax = cloud.Dmax();
        scalar D10 = cloud.Dij(1, 0);
        scalar D32 = cloud.Dij(3, 2);

        Log << type() << " " << name() <<  " write:" << nl
            << "    number of parcels : " << nParcels << nl
            << "    mass in system    : " << massInSystem << nl
            << "    maximum diameter  : " << Dmax << nl
            << "    D10 diameter      : " << D10 << nl
            << "    D32 diameter      : " << D32 << nl
            << endl;

        if (writeToFile())
        {
            writeTime(files(i));
            files(i)
                << token::TAB
                << nParcels << token::TAB
                << massInSystem << token::TAB
                << Dmax << token::TAB
                << D10 << token::TAB
                << D32 << token::TAB
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //
