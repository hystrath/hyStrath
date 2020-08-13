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

#include "LarsenBorgnakkeVT.H"
#include "fvm.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LarsenBorgnakkeVT<ThermoType>::updateCoefficients()
{
    //tauVTijModel_().update();

    const volScalarField& Tt = thermo_.Tt();
    const volScalarField& et = thermo_.het();
    const scalarField& TtCells = Tt.internalField();
    const scalarField& etCells = et.internalField();

    forAll(solvedVibEqSpecies(), speciei)
    {
        /*volScalarField sumMolarWeightedReciprocalRelaxationTimes = 0/tauVTij(0,0);
        volScalarField tmpXe = 0*thermo_.composition().X(0);

        forAll(species(), speciej)
        {
            const volScalarField& Xj = thermo_.composition().X(speciej);

            // no electrons here (Candler 09), although (ScalabrinPhD 07) defined tauVTie
            if (this->thermo_.composition().electronId() != speciej)
            {
                sumMolarWeightedReciprocalRelaxationTimes += Xj/tauVTij(speciei,speciej);
            }
            else
            {
                tmpXe = Xj;
            }
        }

        tauVT_[speciei] = (1.0-tmpXe) / (sumMolarWeightedReciprocalRelaxationTimes
            + dimensionedScalar("SMALL", dimless/dimTime, Foam::SMALL));*/

        const volScalarField& evi = thermo_.composition().hevel(speciei); // shoudl hev but dont know if updated with hyLight
        const volScalarField& zetavi = thermo_.composition().zetav(speciei);

        const scalarField& evCells = evi.internalField();
        const scalarField& zetavCells = zetavi.internalField();

        scalarField& tauVTCells = tauVT_[speciei].internalField();

        scalar thetaV = 3371.0;
        scalar thetaD = 113500.0;
        scalar Tref = 273.0;
        scalar Zref = 52560.0;
        scalar TZref = 3371.0;

        scalar kB = 1.3806e-23;

        scalarField collFrequencyCells(tauVTCells.size(), 0.0);
        scalarField consistencyFactorCells(tauVTCells.size(), 0.0);

        forAll(species(), speciej)
        {
            const volScalarField& nDj = thermo_.composition().nD(speciej);
            const scalarField& nDjCells = nDj.internalField();

            forAll(tauVTCells, celli)
            {
                scalar Tcoll = max(TtCells[celli], thermo_.composition().Tv(speciei)[celli]);

                collFrequencyCells[celli] += nDjCells[celli]*sqr(dref(speciej))
                    *sqrt(8.0*3.1415926*kB*Tref/reducedMolecularMass(speciei, speciej))
                    *pow(Tcoll/Tref, 1.0-omega(speciej));

                consistencyFactorCells[celli] = 1.0 + (0.5*sqr(zetavCells[celli])*exp(thetaV/Tcoll))/(5.0 - 2.0*omega(speciei));
            }

            Info << "collFrequency " << collFrequencyCells << endl;
            Info << "consistencyFactor " << consistencyFactorCells << endl;
        }

        forAll(tauVTCells, celli)
        {
            scalar B = Zref*pow(thetaD/TZref, -omega(speciei));
            scalar D = pow(thetaD/TZref, 0.33333) - 1.0;

            scalar Rs = 2.968e02;

            //scalar Tcoll = floor((2.0*3.0/5.0*etCells[celli] + evCells[celli])/(Rs*thetaV))*thetaV/(3.5-omega); // just transla in etr
            scalar Tcoll = max(TtCells[celli], thermo_.composition().Tv(speciei)[celli]);
            Info << "Ttr " << TtCells[celli] << endl;
            Info << "Tcoll " << Tcoll << endl;

            scalar Zvs = pow(thetaD/Tcoll, omega(speciei))*pow(B, (pow(thetaD/Tcoll, 0.33333) - 1.0)/D);
            Info << "Zvs " << Zvs << endl;

            tauVTCells[celli] = Zvs/collFrequencyCells[celli]; //*consistencyFactorCells[celli];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::LarsenBorgnakkeVT<ThermoType>::LarsenBorgnakkeVT
(
    rho2ReactionThermo& thermo,
    const compressibleTurbulenceModel& turbulence
)
:
    relaxationTimeModel(thermo, turbulence),

    speciesThermo_
    (
        dynamic_cast<const multi2ComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )
{
    tauVT_.setSize(solvedVibEqSpecies().size());

    forAll(tauVT_, speciei)
    {
        tauVT_.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "tauVT_" + solvedVibEqSpecies()[speciei],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("tauVT", dimTime, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::LarsenBorgnakkeVT<ThermoType>::correct()
{
    updateCoefficients();

    const volScalarField& p = thermo_.p();
    const volScalarField& Tt = thermo_.Tt();

    const scalarField& pCells = p.internalField();
    const scalarField& TtCells = Tt.internalField();

    forAll(solvedVibEqSpecies(), speciei)
    {
        const volScalarField& pD = thermo_.composition().pD(speciei);
        const volScalarField& ev = thermo_.composition().hevel(speciei);
        const volScalarField& tauVT = this->tauVT_[speciei];
        volScalarField& QVT = this->QVT_[speciei];


        const scalarField& pDCells = pD.internalField();
        const scalarField& evCells = ev.internalField();
        const scalarField& tauVTCells = tauVT.internalField();
        scalarField& QVTCells = QVT.internalField();

        // Electrons are not included into the calculation. See private member function above.
        forAll(QVTCells, celli)
        {
            const scalar evZCelli = thermo_.composition().HEvel(speciei, pCells[celli], TtCells[celli]);

            QVTCells[celli] = pDCells[celli]/tauVTCells[celli]*(evZCelli - evCells[celli]);
        }

        /*forAll(QVT.boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = p.boundaryField()[patchi];
            const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
            const fvPatchScalarField& ppD = pD.boundaryField()[patchi];
            const fvPatchScalarField& pev = ev.boundaryField()[patchi];
            const fvPatchScalarField& ptauVT = tauVT.boundaryField()[patchi];

            fvPatchScalarField& pQVT = QVT.boundaryField()[patchi];

            forAll(pQVT, facei)
            {
                const scalar pevZFacei = thermo_.composition().HEvel(speciei, pp[facei], pTt[facei]);

                pQVT[facei] = ppD[facei]/ptauVT[facei]*(pevZFacei - pev[facei]);
            }
        }*/
    }
}


template<class ThermoType>
bool Foam::LarsenBorgnakkeVT<ThermoType>::read()
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
