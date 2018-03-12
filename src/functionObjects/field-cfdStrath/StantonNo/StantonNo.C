/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "StantonNo.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(StantonNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        StantonNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::StantonNo::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "St ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::StantonNo::StantonNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    wallHeatFlux_("wallHeatFlux"),
    inflowPatchName_("inlet"),
    wallHeatFluxHeader_
    (
        wallHeatFlux_,
        runTime.timeName(),
        mesh_,
        IOobject::NO_READ
    )
{
    read(dict);

    writeFileHeader(file());

    volScalarField* StantonNoPtr
    (
        new volScalarField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    mesh_.objectRegistry::store(StantonNoPtr);
    
    //- Check if the wallHeatFlux field exists
    wallHeatFluxHeader_ = IOobject
    (
        wallHeatFlux_,
        runTime.timeName(),
        mesh_,
        IOobject::MUST_READ
    );
    
    if (!wallHeatFluxHeader_.typeHeaderOk<volScalarField>(false))
    {
        FatalErrorInFunction
            << "Unable to find the wallHeatFlux"
            << exit(FatalError);
    } 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::StantonNo::~StantonNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::StantonNo::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    
    wallHeatFlux_ = dict.lookupOrDefault<word>("wallHeatFlux", "wallHeatFlux");
        
    inflowPatchName_ = dict.lookupOrDefault<word>("inflowPatchName", "inlet");

    return true;
}


bool Foam::functionObjects::StantonNo::execute()
{
    volScalarField& StantonNo =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(typeName)
        );

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        if (foundObject<multi2Thermo>(multi2Thermo::dictName))
        {    
            volScalarField::Boundary& StantonNoBf =
                StantonNo.boundaryFieldRef();

            const turbulenceModel& model =
                lookupObject<turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );
                
            const multi2Thermo& thermo =
                lookupObject<multi2Thermo>(multi2Thermo::dictName);
                
            const volScalarField wallHeatFlux(wallHeatFluxHeader_, mesh_);
            
            const volScalarField::Boundary& wallHeatFluxBf = 
                wallHeatFlux.boundaryField();

            const label inflowPatchId = 
                mesh_.boundaryMesh().findPatchID(inflowPatchName_);
                
            tmp<volScalarField> trho = thermo.rho();
            const volScalarField::Boundary& rhoBf = trho().boundaryField();

            const volVectorField::Boundary& UBf = model.U().boundaryField();
            
            const scalar rhoinf = rhoBf[inflowPatchId][0];
            const scalar magUinf = mag(UBf[inflowPatchId][0]);

            const fvPatchList& patches = mesh_.boundary();
            
            forAll(patches, patchi)
            {
                const fvPatch& patch = patches[patchi];

                if (isA<wallFvPatch>(patch))
                {
                    StantonNoBf[patchi] =
                        wallHeatFluxBf[patchi]
                       /(0.5*rhoinf*pow(magUinf, 3));
                }
            }
        }
        else
        {
            WarningInFunction
                << "Unable to find thermo model in the "
                << "database: StantonNo will not be calculated" << endl;
            return false;
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: StantonNo will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::StantonNo::write()
{
    const volScalarField& StantonNo =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << StantonNo.name() << endl;

    StantonNo.write();

    const volScalarField::Boundary& StantonNoBf = StantonNo.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& StantonNop = StantonNoBf[patchi];

            const scalar minStatonNo = gMin(StantonNop);
            const scalar maxStatonNo = gMax(StantonNop);
            const scalar avgStatonNo = gAverage(StantonNop);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " St : min = " << minStatonNo << ", max = " << maxStatonNo
                    << ", average = " << avgStatonNo << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minStatonNo
                    << token::TAB << maxStatonNo
                    << token::TAB << avgStatonNo
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
