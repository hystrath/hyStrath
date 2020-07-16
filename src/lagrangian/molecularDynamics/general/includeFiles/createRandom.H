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

//
// createMesh.H
// ~~~~~~~~~~~~
	string inputFile("randomDict_");
	inputFile.append(argv[0]);

	IOdictionary randomDict
    (
        IOobject
        (
			inputFile,
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    //Initialise seed to system time
    label rndSeed = clock::getTime();

    //Determine the size of the random cache by multiplying by the number of molecules
    label rndCacheSizeMultiplier = 20;

    //Read the random seed and over-write system time (if it exists in the dictionary)
    if (randomDict.headerOk())
    {
        if(randomDict.found("randomSeed"))
            rndSeed = readLabel(randomDict.lookup("randomSeed"));

        if(randomDict.found("randomCacheSizeMultiplier"))
            rndCacheSizeMultiplier = readLabel(randomDict.lookup("randomCacheSizeMultiplier"));
    }

    randomDict.close();

    cachedRandomMD rndGen(rndSeed, rndCacheSizeMultiplier);

    IOstream::defaultPrecision(2); //Set write precision for random value output

    fileName outputFile(runTime.path()/"random_"+argv[0]);
    OFstream os(outputFile);

    os << "Random Seed: " << rndSeed << endl;
    os << "Cache Size Multiplier: " << rndCacheSizeMultiplier << endl;

    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    //Read the controlDict file and reset output precision if entry found
    if (controlDict.headerOk())
    {
        if(controlDict.found("writePrecision"))
            IOstream::defaultPrecision(readLabel(controlDict.lookup("writePrecision")));
        else //No default precision specified so reset back to standard default
            IOstream::defaultPrecision(6);
    }

    controlDict.close();
