/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "blendingFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blendingFactor, 0);
    addToRunTimeSelectionTable(functionObject, blendingFactor, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::blendingFactor::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Blending factor");
    writeCommented(os, "Time");
    writeTabbed(os, "Scheme1");
    writeTabbed(os, "Scheme2");
    writeTabbed(os, "Blended");
    os  << endl;
}


bool Foam::functionObjects::blendingFactor::calc()
{
    bool processed = false;

    processed = processed || calcScheme<scalar>();
    processed = processed || calcScheme<vector>();

    return processed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::blendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    phiName_("phi"),
    tolerance_(0.001)
{
    read(dict);
    writeFileHeader(file());
    setResultName(typeName, "");

    tmp<volScalarField> indicatorPtr
    (
        new volScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    store(resultName_, indicatorPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::~blendingFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::blendingFactor::read(const dictionary& dict)
{
    if (fieldExpression::read(dict) && writeFile::read(dict))
    {
        phiName_ = dict.lookupOrDefault<word>("phi", "phi");

        tolerance_ = 0.001;
        if
        (
            dict.readIfPresent("tolerance", tolerance_)
         && (tolerance_ < 0 || tolerance_ > 1)
        )
        {
            FatalErrorInFunction
                << "tolerance must be in the range 0 to 1.  Supplied value: "
                << tolerance_ << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::blendingFactor::write()
{
    if (fieldExpression::write())
    {
        const volScalarField& indicator =
            lookupObject<volScalarField>(resultName_);

        // Generate scheme statistics
        label nCellsScheme1 = 0;
        label nCellsScheme2 = 0;
        label nCellsBlended = 0;
        forAll(indicator, celli)
        {
            scalar i = indicator[celli];

            if (i < tolerance_)
            {
                nCellsScheme1++;
            }
            else if (i > (1 - tolerance_))
            {
                nCellsScheme2++;
            }
            else
            {
                nCellsBlended++;
            }
        }

        reduce(nCellsScheme1, sumOp<label>());
        reduce(nCellsScheme2, sumOp<label>());
        reduce(nCellsBlended, sumOp<label>());

        Log << "    scheme 1 cells :  " << nCellsScheme1 << nl
            << "    scheme 2 cells :  " << nCellsScheme2 << nl
            << "    blended cells  :  " << nCellsBlended << nl
            << endl;

        file()
            << time_.time().value()
            << token::TAB << nCellsScheme1
            << token::TAB << nCellsScheme2
            << token::TAB << nCellsBlended
            << endl;
    }

    return true;
}


// ************************************************************************* //
