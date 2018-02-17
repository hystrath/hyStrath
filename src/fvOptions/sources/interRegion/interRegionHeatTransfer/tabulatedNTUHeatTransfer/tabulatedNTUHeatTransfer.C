/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "tabulatedNTUHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(tabulatedNTUHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            option,
            tabulatedNTUHeatTransfer,
            dictionary
        );
    }

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fv::tabulatedNTUHeatTransfer::geometryModeType,
        2
    >::names[] =
    {
        "calculated",
        "user"
    };
}

const Foam::NamedEnum<Foam::fv::tabulatedNTUHeatTransfer::geometryModeType, 2>
Foam::fv::tabulatedNTUHeatTransfer::geometryModelNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::interpolation2DTable<Foam::scalar>&
Foam::fv::tabulatedNTUHeatTransfer::ntuTable()
{
    if (!ntuTable_.valid())
    {
        ntuTable_.reset(new interpolation2DTable<scalar>(coeffs_));
    }

    return ntuTable_();
}


const Foam::basicThermo& Foam::fv::tabulatedNTUHeatTransfer::thermo
(
    const fvMesh& mesh
) const
{
    if (!mesh.foundObject<basicThermo>("thermophysicalProperties"))
    {
        FatalErrorInFunction
            << " on mesh " << mesh.name()
            << " could not find thermophysicalProperties "
            << exit(FatalError);
    }

    return mesh.lookupObject<basicThermo>("thermophysicalProperties");
}


void Foam::fv::tabulatedNTUHeatTransfer::initialiseGeometry()
{
    if (Ain_ < 0)
    {
        geometryMode_ =
            geometryModelNames_.read(coeffs_.lookup("geometryMode"));

        Info<< "Region " << mesh_.name() << " " << type() << " " << name_ << " "
            << geometryModelNames_[geometryMode_] << " geometry:" << nl;

        switch (geometryMode_)
        {
            case gmCalculated:
            {
                const fvMesh& nbrMesh =
                    mesh_.time().lookupObject<fvMesh>(nbrRegionName());

                word inletPatchName(coeffs_.lookup("inletPatch"));
                word inletPatchNbrName(coeffs_.lookup("inletPatchNbr"));

                Info<< "    Inlet patch           : " << inletPatchName << nl
                    << "    Inlet patch neighbour : " << inletPatchNbrName
                    << nl;

                label patchI = mesh_.boundary().findPatchID(inletPatchName);
                label patchINbr =
                    nbrMesh.boundary().findPatchID(inletPatchNbrName);

                scalar alpha(readScalar(coeffs_.lookup("inletBlockageRatio")));

                if ((alpha < 0) || (alpha > 1))
                {
                    FatalErrorInFunction
                        << "Inlet patch blockage ratio must be between 0 and 1"
                        << ".  Current value: " << alpha
                        << abort(FatalError);
                }

                scalar alphaNbr
                (
                    readScalar(coeffs_.lookup("inletBlockageRatioNbr"))
                );

                if ((alphaNbr < 0) || (alphaNbr > 1))
                {
                    FatalErrorInFunction
                        << "Inlet patch neighbour blockage ratio must be "
                        << "between 0 and 1.  Current value: " << alphaNbr
                        << abort(FatalError);
                }

                Info<< "    Inlet blockage ratio  : " << alpha << nl
                    << "    Inlet blockage ratio neighbour : " << alphaNbr
                    << nl;

                Ain_ =
                    (scalar(1) - alpha)
                   *gSum(mesh_.magSf().boundaryField()[patchI]);

                AinNbr_ =
                    (scalar(1) - alphaNbr)
                   *gSum(nbrMesh.magSf().boundaryField()[patchINbr]);

                scalar beta(readScalar(coeffs_.lookup("coreBlockageRatio")));

                if ((beta < 0) || (beta > 1))
                {
                    FatalErrorInFunction
                        << "Core volume blockage ratio must be between 0 and 1"
                        << ".  Current value: " << beta
                        << abort(FatalError);
                }

                Info<< "    Core volume blockage ratio : " << beta << nl;

                Vcore_ = (scalar(1) - beta)*meshInterp().V();

                break;
            }
            case gmUser:
            {
                coeffs_.lookup("Ain") >> Ain_;
                coeffs_.lookup("AinNbr") >> AinNbr_;

                if (!coeffs_.readIfPresent("Vcore", Vcore_))
                {
                    Vcore_ = meshInterp().V();
                }

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unhandled enumeration " << geometryMode_
                    << abort(FatalError);
            }
        }

        Info<< "    Inlet area local      : " << Ain_ << nl
            << "    Inlet area neighbour  : " << AinNbr_ << nl
            << "    Core volume           : " << Vcore_ << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::tabulatedNTUHeatTransfer::tabulatedNTUHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    UName_(coeffs_.lookupOrDefault<word>("U", "U")),
    UNbrName_(coeffs_.lookupOrDefault<word>("UNbr", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    rhoNbrName_(coeffs_.lookupOrDefault<word>("rhoNbr", "rho")),
    ntuTable_(),
    geometryMode_(gmCalculated),
    Ain_(-1),
    AinNbr_(-1),
    Vcore_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::tabulatedNTUHeatTransfer::~tabulatedNTUHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::tabulatedNTUHeatTransfer::calculateHtc()
{
    initialiseGeometry();

    const fvMesh& nbrMesh = mesh_.time().lookupObject<fvMesh>(nbrRegionName());

    const basicThermo& thermo = this->thermo(mesh_);
    const basicThermo& thermoNbr = this->thermo(nbrMesh);
    const volScalarField Cp(thermo.Cp());
    const volScalarField CpNbr(thermoNbr.Cp());

    // Calculate scaled mass flow for primary region
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName_);
    const scalarField mDot(mag(U)*rho*Ain_);

    // Calculate scaled mass flow for neighbour region
    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);
    const scalarField UMagNbr(mag(UNbr));
    const scalarField UMagNbrMapped(interpolate(UMagNbr));
    const scalarField& rhoNbr =
        nbrMesh.lookupObject<volScalarField>(rhoNbrName_).internalField();
    const scalarField rhoNbrMapped(interpolate(rhoNbr));
    const scalarField mDotNbr(UMagNbrMapped*rhoNbrMapped*AinNbr_);


    scalarField& htcc = htc_.primitiveFieldRef();
    const interpolation2DTable<Foam::scalar>& ntuTable = this->ntuTable();

    forAll(htcc, cellI)
    {
        scalar Cpc = Cp[cellI];
        scalar CpcNbr = CpNbr[cellI];
        scalar mDotc = mDot[cellI];
        scalar mDotcNbr = mDotNbr[cellI];
        scalar Cmin = min(Cpc*mDotc, CpcNbr*mDotcNbr);
        scalar ntu = ntuTable(mDotc, mDotcNbr);

        htcc[cellI] = Cmin*ntu/Vcore_;
    }
}


bool Foam::fv::tabulatedNTUHeatTransfer::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("U", UName_);
        coeffs_.readIfPresent("UNbr", UNbrName_);
        coeffs_.readIfPresent("rho", rhoName_);
        coeffs_.readIfPresent("rhoNbr", rhoNbrName_);

        // Force geometry re-initialisation
        Ain_ = -1;
        initialiseGeometry();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
