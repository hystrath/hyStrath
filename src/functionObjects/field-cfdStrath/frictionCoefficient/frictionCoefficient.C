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

#include "frictionCoefficient.H"
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
    defineTypeNameAndDebug(frictionCoefficient, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        frictionCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::frictionCoefficient::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Cf ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::frictionCoefficient::frictionCoefficient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    wallShearStress_("wallShearStress"),
    inflowPatchName_("inlet"),
    wallShearStressHeader_
    (
        wallShearStress_,
        runTime.timeName(),
        mesh_,
        IOobject::NO_READ
    )
{
    read(dict);

    writeFileHeader(file());

    volScalarField* frictionCoefficientPtr
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

    mesh_.objectRegistry::store(frictionCoefficientPtr);
    
    //- Check if the wallShearStress field exists
    wallShearStressHeader_ = IOobject
    (
        wallShearStress_,
        runTime.timeName(),
        mesh_,
        IOobject::MUST_READ
    );
    
    if (!wallShearStressHeader_.typeHeaderOk<volVectorField>(false))
    {
        FatalErrorInFunction
            << "Unable to find the wallShearStress"
            << exit(FatalError);
    } 
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::frictionCoefficient::~frictionCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::frictionCoefficient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    
    wallShearStress_ = 
        dict.lookupOrDefault<word>("wallShearStress", "wallShearStress");
        
    inflowPatchName_ = dict.lookupOrDefault<word>("inflowPatchName", "inlet");

    return true;
}


bool Foam::functionObjects::frictionCoefficient::execute()
{
    volScalarField& frictionCoefficient =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(typeName)
        );

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        if (foundObject<multi2Thermo>(multi2Thermo::dictName))
        {    
            volScalarField::Boundary& frictionCoefficientBf =
                frictionCoefficient.boundaryFieldRef();

            const turbulenceModel& model =
                lookupObject<turbulenceModel>
                (
                    turbulenceModel::propertiesName
                );
                
            const multi2Thermo& thermo =
                lookupObject<multi2Thermo>(multi2Thermo::dictName);
                
            const volVectorField wallShearStress(wallShearStressHeader_, mesh_);
            
            const volVectorField::Boundary& wallShearStressBf = 
                wallShearStress.boundaryField();

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
                    forAll(patches[patchi], facei)
                    {
                        frictionCoefficientBf[patchi][facei] =
                            mag(wallShearStressBf[patchi][facei])
                           /(0.5*rhoinf*pow(magUinf, 2));
                    }
                }
            }
        }
        else
        {
            WarningInFunction
                << "Unable to find thermo model in the "
                << "database: Cf will not be calculated" << endl;
            return false;
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: Cf will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::frictionCoefficient::write()
{
    const volScalarField& frictionCoefficient =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << frictionCoefficient.name() << endl;

    frictionCoefficient.write();

    const volScalarField::Boundary& frictionCoefficientBf = 
        frictionCoefficient.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& frictionCoefficientp = 
                frictionCoefficientBf[patchi];

            const scalar minCf = gMin(frictionCoefficientp);
            const scalar maxCf = gMax(frictionCoefficientp);
            const scalar avgCf = gAverage(frictionCoefficientp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " Cf : min = " << minCf << ", max = " << maxCf
                    << ", average = " << avgCf << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minCf
                    << token::TAB << maxCf
                    << token::TAB << avgCf
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
