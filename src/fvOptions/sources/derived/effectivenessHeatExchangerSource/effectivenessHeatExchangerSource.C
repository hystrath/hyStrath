/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "effectivenessHeatExchangerSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(effectivenessHeatExchangerSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        effectivenessHeatExchangerSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::initialise()
{
    const label zoneID = mesh_.faceZones().findZoneID(faceZoneName_);

    if (zoneID < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneID];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(facei);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = facei - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::effectivenessHeatExchangerSource::effectivenessHeatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    secondaryMassFlowRate_(0),
    secondaryInletT_(0),
    primaryInletT_(0),
    userPrimaryInletT_(false),
    eTable_(),
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    faceZoneName_("unknown-faceZone"),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict);

    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(1, false);

    eTable_.reset(new interpolation2DTable<scalar>(coeffs_));

    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const surfaceScalarField Tf(fvc::interpolate(T));

    scalar sumPhi = 0;
    scalar sumMagPhi = 0;
    scalar CpfMean = 0;
    scalar primaryInletTfMean = 0;
    forAll(faceId_, i)
    {
        label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            scalar phii = phi.boundaryField()[patchi][facei]*faceSign_[i];

            sumPhi += phii;

            scalar Cpfi = Cpf.boundaryField()[patchi][facei];
            scalar Tfi = Tf.boundaryField()[patchi][facei];
            scalar magPhii = mag(phii);

            sumMagPhi += magPhii;
            CpfMean += Cpfi*magPhii;
            primaryInletTfMean += Tfi*magPhii;
        }
        else
        {
            scalar phii = phi[facei]*faceSign_[i];
            scalar magPhii = mag(phii);

            sumPhi += phii;
            sumMagPhi += magPhii;
            CpfMean += Cpf[facei]*magPhii;
            primaryInletTfMean += Tf[facei]*magPhii;
        }
    }
    reduce(CpfMean, sumOp<scalar>());
    reduce(sumPhi, sumOp<scalar>());
    reduce(sumMagPhi, sumOp<scalar>());
    reduce(primaryInletTfMean, sumOp<scalar>());

    primaryInletTfMean /= sumMagPhi + ROOTVSMALL;
    CpfMean /= sumMagPhi + ROOTVSMALL;

    scalar primaryInletT = primaryInletT_;
    if (!userPrimaryInletT_)
    {
        primaryInletT = primaryInletTfMean;
    }

    scalar Qt =
        eTable_()(mag(sumPhi), secondaryMassFlowRate_)
       *(secondaryInletT_ - primaryInletT)
       *CpfMean*mag(sumPhi);

    const scalarField TCells(T, cells_);
    scalar Tref = 0;
    if (Qt > 0)
    {
        Tref = max(TCells);
        reduce(Tref, maxOp<scalar>());
    }
    else
    {
        Tref = min(TCells);
        reduce(Tref, minOp<scalar>());
    }

    scalarField deltaTCells(cells_.size(), 0);
    forAll(deltaTCells, i)
    {
        if (Qt > 0)
        {
            deltaTCells[i] = max(Tref - TCells[i], 0.0);
        }
        else
        {
            deltaTCells[i] = max(TCells[i] - Tref, 0.0);
        }
    }

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    scalar sumWeight = 0;
    forAll(cells_, i)
    {
        label celli = cells_[i];
        sumWeight += V[celli]*mag(U[celli])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    if (this->V() > VSMALL && mag(Qt) > VSMALL)
    {
        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            heSource[celli] -=
                Qt*V[celli]*mag(U[celli])*deltaTCells[i]
               /(sumWeight + ROOTVSMALL);
        }
    }

    Info<< type() << ": " << name() << nl << incrIndent
        << indent << "Net mass flux [Kg/s]      : " << sumPhi << nl
        << indent << "Total energy exchange [W] : " << Qt << nl
        << indent << "Tref [K]                  : " << Tref << nl
        << indent << "Effectiveness             : "
        << eTable_()(mag(sumPhi), secondaryMassFlowRate_) << decrIndent
        << nl << endl;
}


bool Foam::fv::effectivenessHeatExchangerSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        UName_ = coeffs_.lookupOrDefault<word>("U", "U");
        TName_ = coeffs_.lookupOrDefault<word>("T", "T");
        phiName_ = coeffs_.lookupOrDefault<word>("phi", "phi");
        coeffs_.lookup("faceZone") >> faceZoneName_;

        coeffs_.lookup("secondaryMassFlowRate") >> secondaryMassFlowRate_;
        coeffs_.lookup("secondaryInletT") >> secondaryInletT_;

        if (coeffs_.readIfPresent("primaryInletT", primaryInletT_))
        {
            Info<< type() << " " << this->name() << ": "
                << "employing user-specified primary flow inlet temperature: "
                << primaryInletT_ << endl;

            userPrimaryInletT_ = true;
        }
        else
        {
            Info<< type() << " " << this->name() << ": "
                << "employing flux-weighted primary flow inlet temperature"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
