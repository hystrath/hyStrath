
/**************************************************************************/
/**Setup fields for Electromagnetic solver**/
/*
Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    dimensionedScalar muMag
    (
        transportProperties.lookup("muMag")
    );

    Info<< "Reading field sigma\n" << endl;
    volScalarField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
*/
    Info<< "Reading field phiE\n" << endl;
    volScalarField phiE
    (
        IOobject
        (
            "phiE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
/*
    Info<< "Reading field A\n" << endl;
    volVectorField A
    (
        IOobject
        (
            "A",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Calculating magnetic field B \n" << endl;
    volVectorField B
    (
        IOobject
        (
            "B",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::curl(A)
    );
*/
    Info << "Calculating electric field E \n" << endl;
    volVectorField E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -fvc::grad(phiE) //- fvc::ddt(A)
    );

/*
    Info << "Calculate electrical current density Je \n" << endl;
    volVectorField  Je
    (
        IOobject
        (
            "Je",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -sigma*(fvc::grad(phiE))
    );
*/
/**************************************************************************/
