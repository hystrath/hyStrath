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
    StantonNumber

Description
    Calculates and reports St for all wall patches, for the specified times.

    Default behaviour assumes operating in compressible mode.
    Use the -incompressible option for incompressible RAS cases.
    The wall heat flux must be first outputed into the result directories via
    the utility wallHeatFlux2T

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fluidThermo.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace std;

void calcSt
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& St
)
{
    IOobject wallHeatFluxheader
    (
        "wallHeatFlux",
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
    
    // Check wallHeatFlux and rho exists
    if (wallHeatFluxheader.headerOk() && rhoheader.headerOk())
    {
        volScalarField wallHeatFlux(wallHeatFluxheader, mesh);
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
        
        const volScalarField::GeometricBoundaryField wallPatches =
                    rho.boundaryField();
                    
        forAll(wallPatches, patchi)
        {
          if ( isA<wallFvPatch> ( wallPatches[patchi].patch() ) )
          {     
                
              St = wallHeatFlux / ( 0.5 * rhoInf * (pow(uInfx,3) + pow(uInfy,3) + pow(uInfz,3)) );
              const scalarField& StNumber = St.boundaryField()[patchi];

              Info<< "Patch " << patchi
                  << " named " << wallPatches[patchi].patch().name()
                  << " \e[1;33mSt\e[0m : min: " << gMin(StNumber) << "\e[1;33m max: " << gMax(StNumber)
                  << "\e[0m average: " << gAverage(StNumber) << endl;          
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

        volScalarField St
        (
            IOobject
            (
                "StantonNumber",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("StantonNumber", dimless, 0.0)
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
              
        if (UHeader.headerOk())
        {
            Info<< "Reading field U" << endl;
            volVectorField U(UHeader, mesh);
            calcSt(mesh, runTime, U, St);
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "\e[1;33m" << "Writing St to field " << St.name() << " \e[0m" << endl;

        St.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
