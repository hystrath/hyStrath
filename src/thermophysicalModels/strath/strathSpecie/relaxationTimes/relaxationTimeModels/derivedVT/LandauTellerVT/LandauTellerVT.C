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
    tauVTijModel_().update();
    
    const label electronId = this->thermo_.composition().electronId();

    forAll(solvedVibEqSpecies(), speciei)
    {
        volScalarField sumMolarWeightedReciprocalRelaxationTimes = 0/tauVTij(0,0);
        volScalarField tmpXe = 0.0*thermo_.composition().X(0);

        forAll(species(), speciej)
        {
            const volScalarField& Xj = thermo_.composition().X(speciej);

            // no electrons here (Candler 09), although (ScalabrinPhD 07)
            // defined tauVTie
            if (electronId != speciej)
            {
                sumMolarWeightedReciprocalRelaxationTimes +=
                    Xj/tauVTij(speciei,speciej);
            }
            else
            {
                tmpXe = Xj;
            }
        }

        tauVT_[speciei] =
            (1.0 - tmpXe)
          / (
                sumMolarWeightedReciprocalRelaxationTimes
              + dimensionedScalar("SMALL", dimless/dimTime, Foam::SMALL)
            );
    }


    /*forAll(solvedVibEqSpecies(), i) // TODO ABORTIVE WORK
    {
        forAll(tauVTmode_[i], mi)
        {
          volScalarField tmpRelTime = 0 / tauVTij(0,0);
          volScalarField tmpXe = 0 * thermo_.composition().X(0);
          if(speciesThermo_[i].particleType() > 1)
          {
              forAll(species(), j)
              {
                  const volScalarField& Xj = thermo_.composition().X(j);
                  forAll(tauVTmode_[j], mj)
                  {
                      if(speciesThermo_[j].particleType() > 0)
                      {
                          tmpRelTime += Xj / tauVTij(i,j);
                      }
                      else
                      {
                          tmpXe = Xj;
                      }
                  }
              }
          }
          
          tauVTmode_[i][mi] =
              (1.0-tmpXe)
            / (
                  tmpRelTime 
                + dimensionedScalar("SMALL", dimless/dimTime, Foam::SMALL)
              );
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
    tauVT_.setSize(solvedVibEqSpecies().size());
    //tauVTmode_.setSize(species().size()); // TODO ABORTIVE WORK

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

    /*forAll(tauVTmode_, speciei)  // TODO ABORTIVE WORK
    {
        tauVTmode_.set
        (
            speciei,
            new PtrList<volScalarField>
            (
                thermo_.composition().noVibrationalTemp(speciei)
            )
        );
    }

    forAll(tauVTmode_, speciei)
    {
      forAll(tauVTmode_[speciei], vibMode)
      {
        tauVTmode_[speciei].set
        (
            vibMode,
            new volScalarField
            (
                IOobject
                (
                    "tauVT_" + species()[speciei] + "." + word(vibMode+1),
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
        scalarField& QVTCells = QVT.primitiveFieldRef();

        // Electrons are not included into the calculation.
        // See private member function above.
        forAll(QVTCells, celli)
        {
            const scalar evZCelli =
                thermo_.composition().HEvel
                (
                    speciei,
                    pCells[celli],
                    TtCells[celli]
                );

            QVTCells[celli] =
                pDCells[celli]/tauVTCells[celli]*(evZCelli - evCells[celli]);
        }

        forAll(QVT.boundaryField(), patchi)
        {
            const fvPatchScalarField& pp = p.boundaryField()[patchi];
            const fvPatchScalarField& pTt = Tt.boundaryField()[patchi];
            const fvPatchScalarField& ppD = pD.boundaryField()[patchi];
            const fvPatchScalarField& pev = ev.boundaryField()[patchi];
            const fvPatchScalarField& ptauVT = tauVT.boundaryField()[patchi];

            fvPatchScalarField& pQVT = QVT.boundaryFieldRef()[patchi];

            forAll(pQVT, facei)
            {
                const scalar pevZFacei =
                    thermo_.composition().HEvel(speciei, pp[facei], pTt[facei]);

                pQVT[facei] = ppD[facei]/ptauVT[facei]*(pevZFacei - pev[facei]);
            }
        }

        /*forAll(QVTmode_[speciei], vibMode) // TODO ABORTIVE WORK
        {
            const volScalarField& hvmode =
                thermo_.composition().hevel_mode(speciei, vibMode);
            const scalarField& hvmodeCells = hvmode.internalField();
            const scalarField& tauVTmodeCells =
                this->tauVTmode_[speciei][vibMode].internalField();
            scalarField& QVTmodeCells =
                this->QVTmode_[speciei][vibMode].primitiveFieldRef();

            // electrons are not included into the calculations.
            // See private member function above.
            forAll(QVTmodeCells, celli) 
            {
                scalar hvZmodeCells =
                    thermo_.composition().HEvel_mode
                    (
                        speciei,
                        vibMode,
                        pCells[celli],
                        TtCells[celli]
                    );
                    
                QVTmodeCells[celli] = pDCells[celli]
                    *(hvZmodeCells - hvmodeCells[celli])/tauVTmodeCells[celli];
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
