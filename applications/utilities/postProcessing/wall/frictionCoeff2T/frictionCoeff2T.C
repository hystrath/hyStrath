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
    frictionCoeff2T

Description
    Calculates and reports Cf for all wall patches, for the specified times
    when using RAS2 turbulence models.

    Default behaviour assumes operating in compressible mode.
    Use the -incompressible option for incompressible RAS2 cases.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"

#include "rho2HTCModel.H"
#include "strath/compressible/RAS/RASModel/RAS2Model.H"
#include "mutWallFunction/mutWallFunctionFvPatchScalarField.H"

#include "nearWallDist.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace std;

void calcIncompressibleCf
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& Cf, 
    volScalarField& Cfx, 
    volScalarField& Cfy, 
    volScalarField& Cfz,
    volScalarField& Cn,    
    volScalarField& Cnx, 
    volScalarField& Cny, 
    volScalarField& Cnz, 
    volScalarField& Cw,    
    volScalarField& Cwx, 
    volScalarField& Cwy, 
    volScalarField& Cwz, 
    const bool& laminar
)
{}


void calcCompressibleCf
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& Cf, 
    volScalarField& Cfx, 
    volScalarField& Cfy, 
    volScalarField& Cfz,
    volScalarField& Cn,    
    volScalarField& Cnx, 
    volScalarField& Cny, 
    volScalarField& Cnz, 
    volScalarField& Cw,    
    volScalarField& Cwx, 
    volScalarField& Cwy, 
    volScalarField& Cwz, 
    const bool& laminar
)
{
    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    
    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }
    
    Info<< "Reading field rho" << endl;
    volScalarField rho(rhoHeader, mesh);
    
    IOobject muHeader
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    
    if (!muHeader.headerOk())
    {
        Info<< "    no mu field" << endl;
        return;
    }

    Info<< "Reading field mu" << endl;
    volScalarField mu(muHeader, mesh);

    #include "compressibleCreatePhi.H"

    autoPtr<hTC2Models::rho2HTCModel> reaction
    (
        hTC2Models::rho2HTCModel::New(mesh)
    );
    rho2ReactionThermo& thermo = reaction->thermo();
    
    volScalarField muEff = mu; //thermo.mu();

    if (!laminar)
    {
        autoPtr<compressible::RAS2Model> RASModel
        (
          compressible::RAS2Model::New
          (
              rho,
              U,
              phi,
              thermo
          )
        );
        
        IOobject mutHeader
        (
            "mut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        
        if (!mutHeader.headerOk())
        {
            Info<< "    no mut field" << endl;
            return;
        }

        Info<< "Reading field mut" << endl;
        volScalarField mut(mutHeader, mesh);

        muEff += mut; //RASModel->mut()();     
    }
    
    volSymmTensorField realTau 
    ( 
        -2.0/3.0 * muEff * fvc::div(U)*I //true for mono-atomic gases
        + muEff * twoSymm(fvc::grad(U))  
    );
    
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
    const fvPatchScalarField& fvp3 = rho.boundaryField()[leftI];
    if (fvp3.size())
    {
        rhoLeft = fvp3[0];
    }
    reduce(rhoLeft, maxOp<scalar>());

    dimensionedScalar rhoInf
    (
        "rhoInf",
        dimensionSet(1, -3, 0, 0, 0),
        rhoLeft
    );
             
    const volScalarField::GeometricBoundaryField wallPatches =
        rho.boundaryField();

    forAll(wallPatches, patchi)
    {
      if ( isA<wallFvPatch> ( wallPatches[patchi].patch() ) )
      {  
       
          dimensionedVector myIdentityVector 
          (   
              "myIdentityVector",
              dimless,             
              vector (1.0, 1.0, 1.0)
          );
          dimensionedVector myXVector 
          (   
              "myXVector",
              dimless,             
              vector (1.0, 0.0, 0.0)
          );
          dimensionedVector myYVector 
          (   
              "myYVector",
              dimless,             
              vector (0.0, 1.0, 0.0)
          );
          dimensionedVector myZVector 
          (   
              "myZVector",
              dimless,             
              vector (0.0, 0.0, 1.0)
          );
          
          //Any curved walls   
          surfaceVectorField realTau_w 
          ( 
              -fvc::interpolate( realTau ) & mesh.Sf() / mesh.magSf()
          );
          
          /*surfaceVectorField realTau_w 
          (
            IOobject
            (
                "realTau_w",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("realTau_w", dimensionSet(1, -1, -2, 0 , 0), vector::zero)
          );
          
          realTau_w =
          ( 
              -fvc::interpolate( realTau ) & mesh.Sf() / mesh.magSf()
          );
          
          realTau_w.write();*/
          
          surfaceVectorField normalTau_w 
          ( 
              (realTau_w & mesh.Sf() / mesh.magSf()) * mesh.Sf() / mesh.magSf()
          );
          surfaceVectorField tangentTau_w 
          (
              realTau_w - normalTau_w
          );
             
          
          
          surfaceScalarField fricCoeff 
          ( 
              mag(realTau_w) / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          );
          surfaceScalarField fricNormalCoeff 
          ( 
              mag(normalTau_w) / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricWallCoeff 
          ( 
              mag(tangentTau_w) / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          
          Cf.boundaryField()[patchi] = fricCoeff.boundaryField()[patchi];
          const scalarField& Cfric = Cf.boundaryField()[patchi];
          Cn.boundaryField()[patchi] = fricNormalCoeff.boundaryField()[patchi];
          const scalarField& CNfric = Cn.boundaryField()[patchi];
          Cw.boundaryField()[patchi] = fricWallCoeff.boundaryField()[patchi];
          const scalarField& CWfric = Cw.boundaryField()[patchi];

          Info<< "Patch " << patchi
              << " named " << wallPatches[patchi].patch().name()
              << " \e[1;33mCf\e[0m : min: " << gMin(Cfric) << "\e[1;33m max: " << gMax(Cfric)
              << "\e[0m average: " << gAverage(Cfric) << endl;
          Info<< "Patch " << patchi
              << " named " << wallPatches[patchi].patch().name()
              << " \e[1;33mCn\e[0m : min: " << gMin(CNfric) << "\e[1;33m max: " << gMax(CNfric)
              << "\e[0m average: " << gAverage(CNfric) << endl;
          Info<< "Patch " << patchi
              << " named " << wallPatches[patchi].patch().name()
              << " \e[1;33mCw\e[0m : min: " << gMin(CWfric) << "\e[1;33m max: " << gMax(CWfric)
              << "\e[0m average: " << gAverage(CWfric) << endl;
              
                     
          surfaceScalarField fricCoeffx 
          ( 
              realTau_w & myXVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricCoeffy 
          ( 
              realTau_w & myYVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricCoeffz 
          ( 
              realTau_w & myZVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          );
          surfaceScalarField fricNormalCoeffx 
          ( 
              normalTau_w & myXVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricNormalCoeffy 
          ( 
              normalTau_w & myYVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricNormalCoeffz 
          ( 
              normalTau_w & myZVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          );
          surfaceScalarField fricWallCoeffx 
          ( 
              tangentTau_w & myXVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricWallCoeffy 
          ( 
              tangentTau_w & myYVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          ); 
          surfaceScalarField fricWallCoeffz 
          ( 
              tangentTau_w & myZVector / ( 0.5 * rhoInf * (pow(uInfx,2) + pow(uInfy,2) + pow(uInfz,2)) )
          );
          Cfx.boundaryField()[patchi] = fricCoeffx.boundaryField()[patchi];
          Cfy.boundaryField()[patchi] = fricCoeffy.boundaryField()[patchi];
          Cfz.boundaryField()[patchi] = fricCoeffz.boundaryField()[patchi];
          Cnx.boundaryField()[patchi] = fricNormalCoeffx.boundaryField()[patchi];
          Cny.boundaryField()[patchi] = fricNormalCoeffy.boundaryField()[patchi];
          Cnz.boundaryField()[patchi] = fricNormalCoeffz.boundaryField()[patchi]; 
          Cwx.boundaryField()[patchi] = fricWallCoeffx.boundaryField()[patchi];
          Cwy.boundaryField()[patchi] = fricWallCoeffy.boundaryField()[patchi];
          Cwz.boundaryField()[patchi] = fricWallCoeffz.boundaryField()[patchi];               
             
      }
      else
      {
          Info<< "\e[0;34mPatch " << patchi
          << " named " << wallPatches[patchi].patch().name()
          << " is not a wall\e[0m " << endl;
      } 
    }

}





int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "incompressible",
        "calculate incompressible y+"
    );
    
    argList::addBoolOption
    (
        "laminar",
        "calculate laminar y+"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    const bool incompressible = false; //args.optionFound("incompressible");
    const bool laminar = args.optionFound("laminar");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        
        volScalarField Cf
        (
            IOobject
            (
                "Cf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cf", dimless, 0.0)
        );
        
        volScalarField Cfx
        (
            IOobject
            (
                "Cfx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cfx", dimless, 0.0)
        );
        
        volScalarField Cfy
        (
            IOobject
            (
                "Cfy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cfy", dimless, 0.0)
        );
        
        volScalarField Cfz
        (
            IOobject
            (
                "Cfz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cfz", dimless, 0.0)
        );

        volScalarField Cn
        (
            IOobject
            (
                "Cn",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cn", dimless, 0.0)
        );
        
        volScalarField Cnx
        (
            IOobject
            (
                "Cnx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cnx", dimless, 0.0)
        );
        
        volScalarField Cny
        (
            IOobject
            (
                "Cny",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cny", dimless, 0.0)
        );
        
        volScalarField Cnz
        (
            IOobject
            (
                "Cnz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cnz", dimless, 0.0)
        );

        volScalarField Cw
        (
            IOobject
            (
                "Cw",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cw", dimless, 0.0)
        );
        
        volScalarField Cwx
        (
            IOobject
            (
                "Cwx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cwx", dimless, 0.0)
        );
        
        volScalarField Cwy
        (
            IOobject
            (
                "Cwy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cwy", dimless, 0.0)
        );
        
        volScalarField Cwz
        (
            IOobject
            (
                "Cwz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Cwz", dimless, 0.0)
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        
        char flowRegime[30] = "";
        
        if (UHeader.headerOk())
        {
            Info<< "Reading field U" << endl;
            volVectorField U(UHeader, mesh);


            if (laminar) {strcat(flowRegime,"laminar ");} else {strcat(flowRegime,"turbulent ");}
            if (!incompressible)
            {
                calcCompressibleCf(mesh, runTime, U, Cf, Cfx, Cfy, Cfz, Cn, Cnx, Cny, Cnz, Cw, Cwx, Cwy, Cwz, laminar);
                strcat(flowRegime,"compressible");
            }
            else
            {
                calcIncompressibleCf(mesh, runTime, U, Cf, Cfx, Cfy, Cfz, Cn, Cnx, Cny, Cnz, Cw, Cwx, Cwy, Cwz, laminar);
                strcat(flowRegime,"incompressible");
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }



        Info<< "\e[1;33m" << "Writing Cf to field " << Cf.name() << ", " <<
        flowRegime << " flow \e[0m" << endl;
        Info<< "\e[1;33m" << "Writing Cfx,y,z to fields Cfx,y,z\e[0m" << nl << endl;
        Info<< "\e[1;33m" << "Writing Cn to field " << Cn.name() << ", " <<
        flowRegime << " flow \e[0m" << endl;
        Info<< "\e[1;33m" << "Writing Cnx,y,z to fields Cnx,y,z\e[0m" << nl << endl;
        Info<< "\e[1;33m" << "Writing Cw to field " << Cw.name() << ", " <<
        flowRegime << " flow \e[0m" << endl;
        Info<< "\e[1;33m" << "Writing Cwx,y,z to fields Cwx,y,z\e[0m" << nl << endl;

        Cf.write(); Cfx.write(); Cfy.write(); Cfz.write();
        Cn.write(); Cnx.write(); Cny.write(); Cnz.write();
        Cw.write(); Cwx.write(); Cwy.write(); Cwz.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
