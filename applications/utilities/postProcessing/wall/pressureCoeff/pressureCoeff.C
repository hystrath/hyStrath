/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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

Application
    pressureCoeff

Description
    Calculates and reports Cp for all wall patches, for the specified times.

    Default behaviour assumes operating in compressible mode.
    Use the -incompressible option for incompressible RAS cases.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fluidThermo.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace std;

void calcCp
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& Cp
)
{
    IOobject pheader
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    IOobject rhoheader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    
    // Check p and rho exists
    if (
            pheader.typeHeaderOk<volScalarField>(false)
         && rhoheader.typeHeaderOk<volScalarField>(false)
       )
    {
        volScalarField p(pheader, mesh);
        volScalarField rho(rhoheader, mesh);     
                           
        Info<< "Reading inlet velocity" << endl;
        scalar ULeft = 0.0; scalar VLeft = 0.0; scalar WLeft = 0.0;
        label leftI = mesh.boundaryMesh().findPatchID("inlet");
        const fvPatchVectorField& fvp = U.boundaryField()[leftI];
        if (fvp.size())
        {
            ULeft = fvp[0].x();
            VLeft = fvp[1].y();
            WLeft = fvp[2].z();
        }
        reduce(ULeft, maxOp<scalar>());
        reduce(VLeft, maxOp<scalar>());
        reduce(WLeft, maxOp<scalar>());

        dimensionedScalar uInfx
        (
            "uInfx",
            dimensionSet(0, 1, -1, 0, 0),
            ULeft
        );
        dimensionedScalar uInfy
        (
            "uInfy",
            dimensionSet(0, 1, -1, 0, 0),
            VLeft
        );
        dimensionedScalar uInfz
        (
            "uInfz",
            dimensionSet(0, 1, -1, 0, 0),
            WLeft
        );
          
        Info<< "Components : " << uInfx.value() << " " << uInfy.value() << " " << uInfz.value() << nl << endl;
        
        Info<< "Reading inlet density rhoInf" << endl;
        scalar rhoLeft = 0.0;
        const fvPatchScalarField& fvp2 = rho.boundaryField()[leftI];
        if (fvp2.size())
        {
            rhoLeft = fvp2[0];
        }
        reduce(rhoLeft, maxOp<scalar>());

        dimensionedScalar rhoInf
        (
            "rhoInf",
            dimensionSet(1, -3, 0, 0, 0),
            rhoLeft
        );
        Info<< "rhoInf : " << rhoInf.value() << endl;
        
        Info<< "Reading inlet pressure pInf" << endl;
        scalar pLeft = 0.0;
        const fvPatchScalarField& fvp3 = p.boundaryField()[leftI];
        if (fvp3.size())
        {
            pLeft = fvp3[0];
        }
        reduce(pLeft, maxOp<scalar>());

        dimensionedScalar pInf
        (
            "pInf",
            dimensionSet(1, -1, -2, 0, 0),
            pLeft
        );
        Info<< "pInf : " << pInf.value() << endl;

        const volScalarField::Boundary wallPatches = rho.boundaryField();
                    
        forAll(wallPatches, patchi)
        {
          if ( isA<wallFvPatch> ( wallPatches[patchi].patch() ) )
          {     
                
              Cp = ( p - pInf ) / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) );
              const scalarField& Cpress = Cp.boundaryField()[patchi];

              Info<< "Patch " << patchi
                  << " named " << wallPatches[patchi].patch().name()
                  << " \e[1;33mCp\e[0m : min: " << gMin(Cpress) << "\e[1;33m max: " << gMax(Cpress)
                  << "\e[0m average: " << gAverage(Cpress) << endl;          
          }
          else
          {
              Info<< "\e[0;34mPatch " << patchi
              << " named " << wallPatches[patchi].patch().name()
              << " is not a wall\e[0m " << endl;
          } 
        }
        
    }
    else
    {
        Info<< "    Missing p or rho" << endl;
    }
                  
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"
   
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        volScalarField Cp
        (
            IOobject
            (
                "pressureCoeff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("pressureCoeff", dimless, 0.0)
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
              
        if (UHeader.typeHeaderOk<volScalarField>(false))
        {
            Info<< "Reading field U" << endl;
            volVectorField U(UHeader, mesh);
            calcCp(mesh, runTime, U, Cp);
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "\e[1;33m" << "Writing Cp to field " << Cp.name() << " \e[0m" << endl;

        Cp.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
