/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "multiSpeciesTransportModel.H"
#include "dimensionedConstants.H"
#include "constants.H"

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(multiSpeciesTransportModel, 0);
    defineRunTimeSelectionTable(multiSpeciesTransportModel, fvMesh);
}


// * * * * * * * * * * * * * * Protected Memeber Functions * * * * * * * * * //

void Foam::multiSpeciesTransportModel::calculateJ
(
    const label i
)
{
    if (thermo_.composition().particleType(i) != 0)
    {
        JnonCorrected_[i] = -rhoD(i)*fvc::grad(thermo_.composition().Y(i))
            + JGradp_[i] + JGradT_[i];

        if (solvingForX_)
        {
            const volScalarField Wmix = thermo_.composition().Wmix();

            JnonCorrected_[i] -= rhoD(i)*thermo_.composition().Y(i)
                *fvc::grad(Wmix)/Wmix;
        }
    }
    else
    {
        volVectorField sum = JnonCorrected_[0]
            *thermo_.composition().particleCharge(0)/W(0);

        for(label specier=1 ; specier < species().size(); specier++)
        {
            if (thermo_.composition().particleType(specier) != 0)
            {
                sum +=
                    JnonCorrected_[specier]
                   *thermo_.composition().particleCharge(specier)
                   /W(specier);
            }
        }

        JnonCorrected_[i] = W(i)*sum;
    }
}


void Foam::multiSpeciesTransportModel::calculateSumDiffusionFluxes()
{
    // Uses the non-corrected diffusion fluxes for the calculation
    // of the diffusion fluxes
    sumDiffusionFluxes_ = JnonCorrected_[0];

    for(label speciej=1 ; speciej < species().size(); speciej++)
    {
        if (thermo_.composition().particleType(speciej) != 0)
        {
            sumDiffusionFluxes_ += JnonCorrected_[speciej];
        }
    }
}


Foam::volVectorField
Foam::multiSpeciesTransportModel::Jcorrected(const label i) const
{
    if ((thermo_.composition().particleType(i) != 0) and (not useNonCorrected_))
    {
        return JnonCorrected_[i]
            - thermo_.composition().Y(i)*sumDiffusionFluxes_;
    }
    else
    {
        return JnonCorrected_[i];
    }
}


void Foam::multiSpeciesTransportModel::
pressureGradientContributionToSpeciesMassFlux()
{
    const volVectorField gradLnpToWmix = fvc::grad(thermo_.p())/thermo_.p()
        /thermo_.composition().Wmix();

    forAll(species(), i)
    {
        const dimensionedScalar Wi
        (
            "unitsW",
            dimMass/dimMoles,
            thermo_.composition().W(i)
        );

        JGradp_[i] = -rhoD(i)*gradLnpToWmix*Wi
            *(thermo_.composition().X(i) - thermo_.composition().Y(i));
    }
}


void Foam::multiSpeciesTransportModel::
temperatureGradientContributionToSpeciesMassFlux()
{

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSpeciesTransportModel::multiSpeciesTransportModel
(
    rho2ReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
:
    IOdictionary
    (
        thermo.transportDictionary()
    ),

    mesh_(thermo.Tt().mesh()),
    thermo_(thermo),
    turbulence_(turbulence),

    spMassFlux_(species().size()),
    JnonCorrected_(species().size()),

    sumDiffusionFluxes_
    (
        IOobject
        (
            "sumDiffusionFluxes",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "sumDiffusionFluxes",
            dimMass/dimArea/dimTime, vector::zero
        )
    ),

    JGradp_(species().size()),
    JGradT_(species().size()),

    useNonCorrected_
    (
        subDict("transportModels").subDict("diffusionModelParameters")
            .lookupOrDefault<bool>("useNonCorrectedForm", false)
    ),
    solvingForX_(false),

    addPressureGradientTerm_
    (
        subDict("transportModels").subDict("diffusionModelParameters")
            .lookupOrDefault<bool>("addPressureGradientTerm", false)
    ),
    addTemperatureGradientTerm_
    (
        subDict("transportModels").subDict("diffusionModelParameters")
            .lookupOrDefault<bool>("addTemperatureGradientTerm", false)
    )
{
    const word dictThermoPhy
    (
        fileName(thermo.lookup("foamChemistryThermoFile")).name()
    );

    const word partialModelName =
        word
        (
            thermo.transportDictionary()
                .subDict("transportModels").lookup("multiSpeciesTransport")
        );

    if (partialModelName == "SCEBD")
    {
        solvingForX_ = true;
    }

    DijModel_.set
    (
        new diffusionModel
        (
            IOdictionary::name(),
            dictThermoPhy,
            thermo.p(),
            thermo.pe(),
            thermo.Tt(),
            species()
         )
    );

    forAll(species(), speciei)
    {
        spMassFlux_.set
        (
            speciei,
            new surfaceScalarField
            (
                IOobject
                (
                    "massFlux_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("massFlux", dimMass/dimTime, 0.0)
            )
        );

        JnonCorrected_.set
        (
            speciei,
            new volVectorField
            (
                IOobject
                (
                    "JnonCorrected_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "JnonCorrected",
                    dimMass/dimArea/dimTime,
                    vector::zero
                )
            )
        );

        JGradp_.set
        (
            speciei,
            new volVectorField
            (
                IOobject
                (
                    "JGradp_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "JGradp",
                    dimMass/dimArea/dimTime,
                    vector::zero
                )
            )
        );

        JGradT_.set
        (
            speciei,
            new volVectorField
            (
                IOobject
                (
                    "JGradT_" + species()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "JGradT",
                    dimMass/dimArea/dimTime,
                    vector::zero
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volVectorField
Foam::multiSpeciesTransportModel::multiSpeciesHeatSource() const
{
    volVectorField multiSpeciesHeatSource
    (
        IOobject
        (
            "multiSpeciesHeatSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "multiSpeciesHeatSource",
            dimEnergy/dimArea/dimTime,
            vector::zero
        )
    );

    forAll(species(), speciej)
    {
        if (thermo_.composition().particleType(speciej) != 0)
        {
            const volScalarField pCells = thermo_.p();
            const volScalarField TtCells = thermo_.Tt();
            const volScalarField TvCells = thermo_.composition().Tv(speciej);

            volScalarField haj
            (
                IOobject
                (
                    "haj",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("haj", dimEnergy/dimMass, 0.0)
            );

            forAll(haj, celli)
            {
                haj[celli] =
                    ha
                    (
                        speciej,
                        pCells[celli],
                        TtCells[celli],
                        TvCells[celli]
                    );
            }

            forAll(haj.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp =
                    thermo_.p().boundaryField()[patchi];
                const fvPatchScalarField& pTt =
                    thermo_.Tt().boundaryField()[patchi];
                const fvPatchScalarField& pTv =
                    thermo_.composition().Tv(speciej).boundaryField()[patchi];

                fvPatchScalarField& phaj = haj.boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    phaj[facei] =
                        ha
                        (
                            speciej,
                            pp[facei],
                            pTt[facei],
                            pTv[facei]
                        );
                }
            }

            // The following works whether the corrected or non-corrected form
            // is employed.
            multiSpeciesHeatSource += Jcorrected(speciej)*haj;
        }
    }

    return multiSpeciesHeatSource;
}


void Foam::multiSpeciesTransportModel::getSpeciesMassFlux
(
    const label i,
    const surfaceScalarField& flux
)
{
//    spMassFlux_[i] = flux; // TODO edit this formula
}


Foam::surfaceScalarField
Foam::multiSpeciesTransportModel::getDiffusiveWallHeatFlux() const
{
    const volScalarField& p = thermo_.p();
    const volScalarField& Tt = thermo_.Tt();
    const PtrList<volScalarField>& Tv = thermo_.composition().Tv();
    const PtrList<volScalarField>& Y = thermo_.composition().Y();
    
    surfaceScalarField heatFlux_diff
    (
        IOobject
        (
            "heatFlux_diff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "heatFlux_diff",
            dimensionSet(1, 0, -3, 0, 0),
            0.0
        )
    );
    
    surfaceScalarField sum_diff_fluxes
    (
        IOobject
        (
            "sum_diff_fluxes",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "sum_diff_fluxes",
            dimensionSet(1, -2, -1, 0, 0),
            0.0
        )
    );
    
    PtrList<surfaceScalarField> Is(species().size());
    PtrList<surfaceScalarField> Js(species().size());
                    
    forAll(species(), speciei)
    {
        if (thermo_.composition().particleType(speciei) != 0)
        {
            Is.set
            (
                speciei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "Is_" + Y[speciei].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    -fvc::interpolate(rhoD(speciei))*fvc::snGrad(Y[speciei])
                )
            );

            sum_diff_fluxes += Is[speciei]; 
        }
    }
    
    forAll(species(), speciei)
    {
        if (thermo_.composition().particleType(speciei) != 0)
        {
            Js.set
            (
                speciei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "Is_" + Y[speciei].name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    Is[speciei]
                )
            );
            
            if (not useNonCorrected_)
            {
                Js[speciei] -= fvc::interpolate(Y[speciei])*sum_diff_fluxes;
            }
            
            volScalarField hai
            (
                IOobject
                (
                    "hai",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("hai", dimEnergy/dimMass, 0.0)
            );

            forAll(hai, celli)
            {
                hai[celli] =
                    ha
                    (
                        speciei,
                        p[celli],
                        Tt[celli],
                        Tv[speciei][celli]
                    );
            }

            forAll(hai.boundaryField(), patchi)
            {
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
                const fvPatchScalarField& pTv =
                    thermo_.composition().Tv(speciei).boundaryField()[patchi];

                fvPatchScalarField& phai = hai.boundaryFieldRef()[patchi];

                forAll(pTt, facei)
                {
                    phai[facei] =
                        ha
                        (
                            speciei,
                            pp[facei],
                            pTt[facei],
                            pTv[facei]
                        );
                }
            }  
            
            heatFlux_diff -= Js[speciei]*fvc::interpolate(hai);
        }
    }
    
    return heatFlux_diff;
}


bool Foam::multiSpeciesTransportModel::read()
{
    return regIOobject::read();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
