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

#include "pressureCoefficient.H"
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
    defineTypeNameAndDebug(pressureCoefficient, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        pressureCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void 
Foam::functionObjects::pressureCoefficient::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cp ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressureCoefficient::pressureCoefficient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    inflowPatchName_("inlet")
{
    read(dict);

    writeFileHeader(file());

    volScalarField* pressureCoefficientPtr
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

    mesh_.objectRegistry::store(pressureCoefficientPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pressureCoefficient::~pressureCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressureCoefficient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    
    inflowPatchName_ = dict.lookupOrDefault<word>("inflowPatchName", "inlet");

    return true;
}


bool Foam::functionObjects::pressureCoefficient::execute()
{
    volScalarField& pressureCoefficient =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(typeName)
        );

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        if (foundObject<multi2Thermo>(multi2Thermo::dictName))
        {    
            volScalarField::Boundary& pressureCoefficientBf =
                pressureCoefficient.boundaryFieldRef();

            const turbulenceModel& model =
                lookupObject<turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );
                
            const multi2Thermo& thermo =
                lookupObject<multi2Thermo>(multi2Thermo::dictName);
                
            const volScalarField::Boundary& pBf = thermo.p().boundaryField();

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
                    pressureCoefficientBf[patchi] =
                        pBf[patchi]/(0.5*rhoinf*pow(magUinf, 2));
                }
            }
        }
        else
        {
            WarningInFunction
                << "Unable to find thermo model in the "
                << "database: Cp will not be calculated" << endl;
            return false;
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: Cp will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::pressureCoefficient::write()
{
    const volScalarField& pressureCoefficient =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << pressureCoefficient.name() << endl;

    pressureCoefficient.write();

    const volScalarField::Boundary& pressureCoefficientBf = 
        pressureCoefficient.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& pressureCoefficientp = 
                pressureCoefficientBf[patchi];

            const scalar minCp = gMin(pressureCoefficientp);
            const scalar maxCp = gMax(pressureCoefficientp);
            const scalar avgCp = gAverage(pressureCoefficientp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " Cp : min = " << minCp << ", max = " << maxCp
                    << ", average = " << avgCp << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minCp
                    << token::TAB << maxCp
                    << token::TAB << avgCp
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
