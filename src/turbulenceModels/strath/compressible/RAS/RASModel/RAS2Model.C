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

\*---------------------------------------------------------------------------*/

#include "RAS2Model.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RAS2Model, 0);
defineRunTimeSelectionTable(RAS2Model, dictionary);
addToRunTimeSelectionTable(turbulenceModel2, RAS2Model, turbulenceModel2);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RAS2Model::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RAS2Model::RAS2Model
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const rho2ReactionThermo& thermophysicalModel, // NEW VINCENT 11/03/2016
    const word& turbulenceModel2Name
)
:
    turbulenceModel2(rho, U, phi, thermophysicalModel, turbulenceModel2Name),

    IOdictionary
    (
        IOobject
        (
            "RASProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    turbulence_(lookup("turbulence")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),

    kMin_("kMin", sqr(dimVelocity), SMALL),
    epsilonMin_("epsilonMin", kMin_.dimensions()/dimTime, SMALL),
    omegaMin_("omegaMin", dimless/dimTime, SMALL)
{
    kMin_.readIfPresent(*this);
    epsilonMin_.readIfPresent(*this);
    omegaMin_.readIfPresent(*this);

    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<RAS2Model> RAS2Model::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const rho2ReactionThermo& thermophysicalModel, // NEW VINCENT 11/03/2016
    const word& turbulenceModel2Name
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "RASProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("RAS2Model")
    );

    Info<< "Selecting RAS turbulence model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RAS2Model::New"
            "("
                "const volScalarField&, "
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "rho2ReactionThermo&, "
                "const word&"
            ")"
        )   << "Unknown RAS2Model type "
            << modelType << nl << nl
            << "Valid RAS2Model types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RAS2Model>
    (
        cstrIter()(rho, U, phi, thermophysicalModel, turbulenceModel2Name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void RAS2Model::correct()
{
    turbulenceModel2::correct();
}


bool RAS2Model::read()
{
    //if (regIOobject::read())

    // Bit of trickery : we are both IOdictionary ('RASProperties') and
    // an regIOobject from the turbulenceModel2 level. Problem is to distinguish
    // between the two - we only want to reread the IOdictionary.

    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            IOdictionary::type()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        lookup("turbulence") >> turbulence_;

        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffDict_ <<= *dictPtr;
        }

        kMin_.readIfPresent(*this);
        epsilonMin_.readIfPresent(*this);
        omegaMin_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
