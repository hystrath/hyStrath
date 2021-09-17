/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

    dictionary squareZoneDict =
        IOdictionary
        (
            IOobject
            (
                "zonesDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

PtrList<entry> regionList(squareZoneDict.lookup("cellZones"));
List<dictionary> dictionaries(regionList.size());
List<word> regionZoneNames(regionList.size());

forAll(regionList, sL)
{
    const entry& region = regionList[sL];
    const dictionary& regionDict = region.dict();
    dictionaries[sL] = regionDict;
    
    const word zoneName(regionDict.lookup("zoneName"));

    if (findIndex(regionZoneNames, zoneName) != -1)
    {
        FatalErrorIn(args.executable())
            << "Different zone names must be used: " 
            << zoneName << nl << "in: "
            << args.caseName()/runTime.system()/"zonesDict"
            << exit(FatalError);
    }
    
    regionZoneNames[sL] = zoneName;
}

