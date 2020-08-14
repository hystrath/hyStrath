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

#include "heRho2Thermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsi2Thermo, class MixtureType>
void Foam::heRho2Thermo<BasicPsi2Thermo, MixtureType>::init2()
{
    //- This function is executed before the initialisation in
    // rho2ReactionThermo. The minimum is done here to calculate
    // the chemical fractions.
    const scalarField& pCells = this->p_.internalField();
    const scalarField& TtCells = this->T_.internalField();

    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();

    forAll(TtCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        psiCells[celli] = mixture_.psi(pCells[celli], TtCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TtCells[celli]);
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pTt = this->T_.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppsi = this->psi_.boundaryFieldRef()[patchi];
        fvPatchScalarField& prho = this->rho_.boundaryFieldRef()[patchi];

        fvPatchScalarField& pht = this->het_.boundaryFieldRef()[patchi];

        if (pTt.fixesValue())
        {
            forAll(pTt, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                ppsi[facei] = mixture_.psi(pp[facei], pTt[facei]);
                prho[facei] = mixture_.rho(pp[facei], pTt[facei]);
            }
        }
        else
        {
            forAll(pTt, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pTt[facei] = mixture_.TtHEt(pht[facei], pp[facei], pTt[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pTt[facei]);
                prho[facei] = mixture_.rho(pp[facei], pTt[facei]);
            }
        }
    }
}


template<class BasicPsi2Thermo, class MixtureType>
void Foam::heRho2Thermo<BasicPsi2Thermo, MixtureType>::calculate2()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsi2Thermo, class MixtureType>
Foam::heRho2Thermo<BasicPsi2Thermo, MixtureType>::heRho2Thermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    he2Thermo<BasicPsi2Thermo, MixtureType>(mesh, phaseName)
{
    init2();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsi2Thermo, class MixtureType>
Foam::heRho2Thermo<BasicPsi2Thermo, MixtureType>::~heRho2Thermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsi2Thermo, class MixtureType>
void Foam::heRho2Thermo<BasicPsi2Thermo, MixtureType>::correct()
{
    if (debug)
    {
        Info<< "entering heRho2Thermo<MixtureType>::correct()" << endl;
    }

    //calculate2();

    if (debug)
    {
        Info<< "exiting heRho2Thermo<MixtureType>::correct()" << endl;
    }
}


// ************************************************************************* //
