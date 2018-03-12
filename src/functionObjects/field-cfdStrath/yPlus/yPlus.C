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

#include "yPlus.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(yPlus, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        yPlus,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::yPlus::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+ ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::yPlus
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict)
{
    read(dict);

    writeFileHeader(file());

    volScalarField* yPlusPtr
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

    mesh_.objectRegistry::store(yPlusPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::yPlus::~yPlus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::yPlus::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    return true;
}


bool Foam::functionObjects::yPlus::execute()
{
    volScalarField& yPlus =
        const_cast<volScalarField&>
        (
            lookupObject<volScalarField>(typeName)
        );

    if (foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        volScalarField::Boundary& yPlusBf = yPlus.boundaryFieldRef();

        const turbulenceModel& model =
            lookupObject<turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        const nearWallDist nwd(mesh_);
        const volScalarField::Boundary& d = nwd.y();

        // nut needed for wall function patches
        tmp<volScalarField> tnut = model.nut();
        const volScalarField::Boundary& nutBf = tnut().boundaryField();

        // U needed for plain wall patches
        const volVectorField::Boundary& UBf = model.U().boundaryField();

        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
            const fvPatch& patch = patches[patchi];

            if (isA<nutWallFunctionFvPatchScalarField>(nutBf[patchi]))
            {
                const nutWallFunctionFvPatchScalarField& nutPf =
                    dynamic_cast<const nutWallFunctionFvPatchScalarField&>
                    (
                        nutBf[patchi]
                    );

                yPlusBf[patchi] = nutPf.yPlus();
            }
            else if (isA<wallFvPatch>(patch))
            {
                yPlusBf[patchi] =
                    d[patchi]
                   *sqrt
                    (
                        model.nuEff(patchi)
                       *mag(UBf[patchi].snGrad())
                    )/model.nu(patchi);
            }
        }
    }
    else
    {
        WarningInFunction
            << "Unable to find turbulence model in the "
            << "database: yPlus will not be calculated" << endl;
        return false;
    }

    return true;
}


bool Foam::functionObjects::yPlus::write()
{
    const volScalarField& yPlus =
        obr_.lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << yPlus.name() << endl;

    yPlus.write();

    const volScalarField::Boundary& yPlusBf = yPlus.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& yPlusp = yPlusBf[patchi];

            const scalar minYplus = gMin(yPlusp);
            const scalar maxYplus = gMax(yPlusp);
            const scalar avgYplus = gAverage(yPlusp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " y+ : min = " << minYplus << ", max = " << maxYplus
                    << ", average = " << avgYplus << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minYplus
                    << token::TAB << maxYplus
                    << token::TAB << avgYplus
                    << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
