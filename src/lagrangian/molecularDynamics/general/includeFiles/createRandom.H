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
