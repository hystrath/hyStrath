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

#include "KnabVV.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::KnabVV<ThermoType>::updateCoefficients()
{
    tauVVijModel_().update();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::KnabVV<ThermoType>::KnabVV
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    relaxationTimeModelVV(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{
    tauVV_.setSize(solvedVibEqSpecies().size());

    forAll(tauVV_, speciei)
    {
        tauVV_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "tauVV_" + solvedVibEqSpecies()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("tauVV", dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::KnabVV<ThermoType>::correct()
{
    updateCoefficients();

    forAll(solvedVibEqSpecies(), speciei)
    {
        const volScalarField& Tt = thermo_.T();
        const volScalarField& p = thermo_.p();
        const volScalarField& pD = thermo_.composition().pD(speciei);
        const volScalarField& ev = thermo_.composition().hevel(speciei);
        volScalarField& QVV = this->QVV_[speciei];

        const scalarField& pCells = p.internalField();
        const scalarField& TtCells = Tt.internalField();
        const scalarField& evCells = ev.internalField();
        scalarField& QVVCells = QVV.primitiveFieldRef();

        QVVCells = 0.0;

        forAll(solvedVibEqSpecies(), speciej)
        {
            if(speciei != speciej)
            {
                const volScalarField& pDj = thermo_.composition().pD(speciej);
                const volScalarField& evj = thermo_.composition().hevel(speciej);
                const volScalarField& factorVVij = tauVVij(speciei, speciej);

                const scalarField& pDjCells = pDj.internalField();
                const scalarField& evjCells = evj.internalField();
                const scalarField& factorVVijCells = factorVVij.internalField();

                forAll(QVVCells, celli)
                {
                    const scalar evZiCelli = thermo_.composition().HEvel(speciei, pCells[celli], TtCells[celli]);
                    const scalar evZjCelli = thermo_.composition().HEvel(speciej, pCells[celli], TtCells[celli]);

                    QVVCells[celli] += factorVVijCells[celli]*pDjCells[celli]/W(speciej)
                        *(evZiCelli/evZjCelli*evjCells[celli] - evCells[celli]);
                }

                forAll(QVV.boundaryField(), patchi)
                {
                    const fvPatchScalarField& pp = p.boundaryField()[patchi];
                    const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
                    const fvPatchScalarField& ppDj = pDj.boundaryField()[patchi];
                    const fvPatchScalarField& pev = ev.boundaryField()[patchi];
                    const fvPatchScalarField& pevj = evj.boundaryField()[patchi];
                    const fvPatchScalarField& pfactorVVij = factorVVij.boundaryField()[patchi];

                    fvPatchScalarField& pQVV = QVV.boundaryFieldRef()[patchi];
                    pQVV = 0;

                    forAll(pQVV, facei)
                    {
                        const scalar pevZiFacei = thermo_.composition().HEvel(speciei, pp[facei], pTt[facei]);
                        const scalar pevZjFacei = thermo_.composition().HEvel(speciej, pp[facei], pTt[facei]);

                        pQVV[facei] += pfactorVVij[facei]*ppDj[facei]/W(speciej)
                            *(pevZiFacei/pevZjFacei*pevj[facei] - pev[facei]);
                    }
                }
            }
        }

        QVV.primitiveFieldRef() *= constant::physicoChemical::NA.value()*pD.internalField()
            *sqrt(8.0*constant::physicoChemical::R.value()*Tt.internalField()/constant::mathematical::pi);

        QVV.boundaryFieldRef() *= constant::physicoChemical::NA.value()*pD.boundaryField()
            *sqrt(8.0*constant::physicoChemical::R.value()*Tt.boundaryField()/constant::mathematical::pi);
    }
}


template<class ThermoType>
bool Foam::KnabVV<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
