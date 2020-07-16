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

\*---------------------------------------------------------------------------*/

reducedUnits refUnits;

IOobject reducedUnitsDictIOobject
(
    "reducedUnitsDict",
    runTime.system(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

if (reducedUnitsDictIOobject.headerOk())
{
    Info<< nl
        << "Reading reference quantities from reducedUnitsDict file." << endl;

    IOdictionary reducedUnitsDict(reducedUnitsDictIOobject);

    refUnits.setRefValues(reducedUnitsDict);
}

Info << refUnits << endl;
