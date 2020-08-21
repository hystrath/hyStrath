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

#include "LandauTellerVT.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::LandauTellerVT<ThermoType>::updateCoefficients()
{
    const dimensionedScalar zero =
        dimensionedScalar("zero", dimless/dimTime, 0.0);
        
    const dimensionedScalar invtSMALL =
        dimensionedScalar("SMALL", dimless/dimTime, Foam::SMALL);
    
    const label electronId = this->thermo_.composition().electronId();
    
    tmp<volScalarField> tXe
    (
        new volScalarField
        (
            IOobject
            (
                "Xe",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    if (electronId != -1)
    {
        volScalarField& XeCells = tXe.ref();
        XeCells = thermo_.composition().X(electronId);
    }
    
    tmp<volScalarField> tSumWeightedReciprocalRelaxTimes
    (
        new volScalarField
        (
            IOobject
            (
                "sumWeightedReciprocalRelaxTimes",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimTime, 0.0)
        )
    );
    
    volScalarField& sumWeightedReciprocalRelaxTimes =
        tSumWeightedReciprocalRelaxTimes.ref();
    
    tauVTijModel_().update();
    
    forAll(molecules(), moli)
    {
        sumWeightedReciprocalRelaxTimes = zero;
        
        // no electrons here (Candler 09), although (ScalabrinPhD 07)
        // defined tauVTie
        forAll(thermo_.composition().heavySpecies(), j)
        {
            const label speciej = thermo_.composition().heavySpeciesIds(j);
            
            const volScalarField& Xj = thermo_.composition().X(speciej);

            sumWeightedReciprocalRelaxTimes += Xj / tauVTij(moli,speciej);
        }

        tauVT_[moli] = (1.0 - tXe.ref())
            / (sumWeightedReciprocalRelaxTimes + invtSMALL);
    }


    /*forAll(molecules(), moli) // ABORTIVE WORK
    {
        forAll(tauVTmode_[moli], mi)
        {
            sumWeightedReciprocalRelaxTimes = zero;
            
            forAll(thermo_.composition().heavySpecies(), j)
            {
                const label speciej = thermo_.composition().heavySpeciesIds(j);
                
                const volScalarField& Xj = thermo_.composition().X(speciej);
                
                forAll(tauVTmode_[speciej], mj)
                {
                    sumWeightedReciprocalRelaxTimes +=
                        Xj / tauVTij(moli,speciej); // EDIT
                }
            }
            
            tauVTmode_[moli][mi] = (1.0 - tXe.ref())
                / (sumWeightedReciprocalRelaxTimes + invtSMALL);
        }
    }*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::LandauTellerVT<ThermoType>::LandauTellerVT
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
    tauVT_.setSize(molecules().size());
    //tauVTmode_.setSize(molecules().size()); // ABORTIVE WORK

    forAll(tauVT_, moli)
    {
        tauVT_.set
        (
            moli,
            new volScalarField
            (
                IOobject
                (
                    "tauVT_" + molecules()[moli],
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

    /*forAll(tauVTmode_, moli)  // ABORTIVE WORK
    {
        tauVTmode_.set
        (
            moli,
            new PtrList<volScalarField>
            (
                thermo_.composition().noVibrationalTemp(moli)
            )
        );
    }

    forAll(tauVTmode_, moli)
    {
        forAll(tauVTmode_[moli], vibMode)
        {
            tauVTmode_[moli].set
            (
                vibMode,
                new volScalarField
                (
                    IOobject
                    (
                        "tauVT_" + molecules()[moli] + "." + word(vibMode+1),
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
    }*/
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::LandauTellerVT<ThermoType>::correct()
{
    updateCoefficients();

    const volScalarField& p = thermo_.p();
    const volScalarField& T = thermo_.T();

    const scalarField& pCells = p.internalField();
    const scalarField& TCells = T.internalField();

    forAll(molecules(), moli)
    {
        const label speciei = thermo_.composition().moleculeIds(moli);
        
        const volScalarField& pD = thermo_.composition().pD(speciei);
        const volScalarField& ev = thermo_.composition().hevel(speciei);
        const volScalarField& tauVT = this->tauVT_[moli];
        volScalarField& QVT = this->QVT_[moli];

        const scalarField& pDCells = pD.internalField();
        const scalarField& evCells = ev.internalField();
        const scalarField& tauVTCells = tauVT.internalField();
        scalarField& QVTCells = QVT.primitiveFieldRef();

        forAll(QVTCells, celli)
        {
            const scalar evZCelli =
                thermo_.composition().HEvel
                (
                    speciei,
                    pCells[celli],
                    TCells[celli]
                );

            QVTCells[celli] = pDCells[celli] / tauVTCells[celli]
                *(evZCelli - evCells[celli]);
        }

        /*forAll(QVTmode_[moli], vibMode) // ABORTIVE WORK
        {
            const volScalarField& hvmode =
                thermo_.composition().hevel_mode(speciei, vibMode);
            const scalarField& hvmodeCells = hvmode.internalField();
            const scalarField& tauVTmodeCells =
                this->tauVTmode_[moli][vibMode].internalField();
            scalarField& QVTmodeCells =
                this->QVTmode_[moli][vibMode].primitiveFieldRef();

            forAll(QVTmodeCells, celli) 
            {
                scalar hvZmodeCells =
                    thermo_.composition().HEvel_mode
                    (
                        speciei,
                        vibMode,
                        pCells[celli],
                        TCells[celli]
                    );
                    
                QVTmodeCells[celli] = pDCells[celli] / tauVTmodeCells[celli]
                    *(hvZmodeCells - hvmodeCells[celli]);
            }
        }*/
    }
}


template<class ThermoType>
bool Foam::LandauTellerVT<ThermoType>::read()
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
