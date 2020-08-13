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

#include "basic2MultiComponentMixture.H"

#include "thermoPhysics2Types.H"

#include "zeroGradientFvPatchFields.H"

#include "fixed2VELEnergyFvPatchScalarField.H"
#include "gradient2VELEnergyFvPatchScalarField.H"
#include "mixed2VELEnergyFvPatchScalarField.H"
#include "fixed2EnergyFvPatchScalarField.H"
#include "gradient2EnergyFvPatchScalarField.H"
#include "mixed2EnergyFvPatchScalarField.H"

#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"

#include "IFstream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basic2MultiComponentMixture::hev2BoundaryBaseTypes
(
    const label speciei,
    const volScalarField::Boundary& tbf1T,
    const bool downgradeToSingleTemperature
)
{
    volScalarField::Boundary tbf =
        this->spTv_[speciei].boundaryField();

    if (downgradeToSingleTemperature)
    {
        tbf = tbf1T;
    }

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList Foam::basic2MultiComponentMixture::he2BoundaryBaseTypes
(
    const volScalarField::Boundary& tbf
)
{
    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
        if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);

            hbt[patchi] = pf.interfaceFieldType();
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            const fixedJumpAMIFvPatchScalarField& pf =
                dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
                (
                    tbf[patchi]
                );

            hbt[patchi] = pf.interfaceFieldType();
        }
    }

    return hbt;
}


Foam::wordList Foam::basic2MultiComponentMixture::hev2BoundaryTypes
(
    const label speciei,
    const volScalarField::Boundary& tbf1T,
    const bool downgradeToSingleTemperature
)
{
    volScalarField::Boundary tbf =
        this->spTv_[speciei].boundaryField();

    if (downgradeToSingleTemperature)
    {
        tbf = tbf1T;
    }

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixed2VELEnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradient2VELEnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixed2VELEnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
        else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
        {
            hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
        }
    }

    return hbt;
}


Foam::wordList Foam::basic2MultiComponentMixture::he2BoundaryTypes
(
    const volScalarField::Boundary& tbf,
    const bool downgradeToSingleTemperature
)
{
    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixed2EnergyFvPatchScalarField::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            //hbt[patchi] = fixed2EnergyFvPatchScalarField::typeName;
            hbt[patchi] = gradient2EnergyFvPatchScalarField::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            /*if (downgradeToSingleTemperature)
            {
                hbt[patchi] = mixed2EnergyFvPatchScalarField::typeName;
            }
            else
            {*/
                hbt[patchi] = fixed2EnergyFvPatchScalarField::typeName;
            //}
            //hbt[patchi] = mixed2EnergyFvPatchScalarField::typeName;
        }
        else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpFvPatchScalarField::typeName;

        }
        else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = energyJumpAMIFvPatchScalarField::typeName;
        }
        else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
        {
            hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
        }
    }

    return hbt;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basic2MultiComponentMixture::basic2MultiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh
)
:
    species_(specieNames),
    molecules_(),
    moleculeIds_(),
    atoms_(),
    atomIds_(),
    ions_(),
    ionIds_(),
    heavySpecies_(),
    heavySpeciesIds_(),
    electronId_(-1),
    active_(species_.size(), true),
    solvedVibEqSpecies_(),
    vibTempAssociativity_(species_.size()),
    Y_(species_.size()),
    X_(species_.size()),
    nD_(species_.size()),
    pP_(species_.size()),
    pD_(species_.size()),
    spTv_(species_.size()),
    //spmodeTv_(species_.size()), TODO ABORTIVE WORK
    hev_(species_.size()),
    heel_(species_.size()),
    hevel_(species_.size()),
    h_(species_.size()),
    //modehevel_(species_.size()), TODO ABORTIVE WORK
    zetaRot_(species_.size()),
    zetaVib_(species_.size()),
    //modezetaVib_(species_.size()), TODO ABORTIVE WORK
    zetaElec_(species_.size()),
    Wmix_
    (
        IOobject
        (
            "Wmix",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Wmix", dimMass/dimMoles, 0.0)
    )
{
    const dictionary thermoDEM =
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
        )()
    );
    
    IOobject TtHeader
    (
        "Tt",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    volScalarField Tt(TtHeader, mesh);
    bool downgradeToSingleTemperature = false;
    bool downgradeToSingleTv = true;

    forAll(species_, i)
    {
        IOobject YHeader
        (
            species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        IOobject XHeader
        (
            "X_" + species_[i],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        );

        bool headersOK = false;

        // check if field exists and can be read
        if
        (
            YHeader.typeHeaderOk<volScalarField>(false)
         && TtHeader.typeHeaderOk<volScalarField>(false)
        )
        {
            headersOK = true;

            vibTempAssociativity_.set
            (
                i,
                new label(0)
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh
                )
            );

            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("X_" + species_[i], dimless, 0.0)
                )
            );
        }
        else if
        (
            XHeader.typeHeaderOk<volScalarField>(false)
         && TtHeader.typeHeaderOk<volScalarField>(false)
        )
        {
            headersOK = true;

            X_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "X_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );

            Y_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(species_[i], dimless, 0.0)
                )
            );
        }

        if (headersOK)
        {
            nD_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "nD_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "nD_" + species_[i],
                        dimless/dimVolume,
                        0.0
                    )
                )
            );

            pP_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "p_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("p_" + species_[i], dimPressure, 0.0)
                )
            );

            pD_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "rho_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("rho_" + species_[i], dimDensity, 0.0)
                )
            );

            {
                IOdictionary thermophysicalProperties
                (
                    IOobject
                    (
                        "thermophysicalProperties",
                        mesh.time().constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                downgradeToSingleTemperature =
                    thermophysicalProperties.lookupOrDefault<bool>
                    (
                        "downgradeToSingleTemperature",
                        false
                    );
                downgradeToSingleTv =
                    thermophysicalProperties.lookupOrDefault<bool>
                    (
                        "downgradeToSingleTv",
                        true
                    );
            }

            if (downgradeToSingleTemperature or downgradeToSingleTv)
            {
                spTv_.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Tv_" + species_[i],
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimTemperature
                    )
                );

                spTv_[i] = Tt;
            }
            else
            {
                spTv_.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Tv_" + species_[i],
                            mesh.time().timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh
                    )
                );
            }

            // TODO ABORTIVE WORK
            /*const label noVibModes =
                readScalar
                (
                    thermoDEM.subDict(species_[i]).subDict("specie").lookup
                    (
                        "noVibTemp"
                    )
                );

            PtrList<volScalarField> fillList1(noVibModes);
            PtrList<volScalarField> fillList2(noVibModes);
            PtrList<volScalarField> fillList3(noVibModes);

            for(label vibMode = 0 ; vibMode < noVibModes ; vibMode++)
            {
                fillList1.set
                (
                    vibMode,
                    new volScalarField
                    (
                        IOobject
                        (
                            "Tv_" + species_[i] + "." + name(vibMode+1),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "Tv_" + species_[i] + "." + name(vibMode+1),
                            dimTemperature,
                            0.0
                        )
                    )
                );

                fillList2.set
                (
                    vibMode,
                    new volScalarField
                    (
                        IOobject
                        (
                            "hevel_" + species_[i] + "." + name(vibMode+1),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "hevel_" + species_[i] + "." + name(vibMode+1),
                            dimEnergy/dimMass,
                            0.0
                        ),
                        this->hev2BoundaryTypes(i, Tt.boundaryField()), // TODO wrong for now introduce mode in bdry cdt
                        this->hev2BoundaryBaseTypes(i, Tt.boundaryField()) // TODO wrong for now
                    )
                );

                fillList3.set
                (
                    vibMode,
                    new volScalarField
                    (
                        IOobject
                        (
                            "zetav_" + species_[i] + "." + name(vibMode+1),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "zetav_" + species_[i] + "." + name(vibMode+1),
                            dimless,
                            0.0
                        )
                    )
                );
            }

            spmodeTv_.set
            (
                i, new PtrList<volScalarField>(noVibModes)
            );

            modehevel_.set
            (
                i, new PtrList<volScalarField>(noVibModes)
            );

            modezetaVib_.set
            (
                i, new PtrList<volScalarField>(noVibModes)
            );

            forAll(spmodeTv_[i], mode)
            {
                spmodeTv_[i].set
                (
                    mode, fillList1[mode]
                );

                modehevel_[i].set
                (
                    mode, fillList2[mode]
                );

                modezetaVib_[i].set
                (
                    mode, fillList3[mode]
                );
            }*/
            // END ABORTIVE WORK

            hev_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "ev_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "ev_" + species_[i],
                        dimEnergy/dimMass,
                        0.0
                    )
                )
            );

            heel_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "eel_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "eel_" + species_[i],
                        dimEnergy/dimMass,
                        0.0
                    )
                )
            );

            if (downgradeToSingleTemperature or downgradeToSingleTv)
            {
                hevel_.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject
                        (
                            "evel_" + species_[i],
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "evel_" + species_[i],
                            dimEnergy/dimMass,
                            0.0
                        )
                    )
                );
            }
            else
            {
                hevel_.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject
                        (
                            "evel_" + species_[i],
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "evel_" + species_[i],
                            dimEnergy/dimMass,
                            0.0
                        ),
                        this->hev2BoundaryTypes(i, Tt.boundaryField()),
                        this->hev2BoundaryBaseTypes(i, Tt.boundaryField())
                    )
                );
            }

            h_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "h_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar
                    (
                        "h_" + species_[i],
                        dimEnergy/dimMass,
                        0.0
                    )
                )
            );

            zetaRot_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "zetar_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zetar_" + species_[i], dimless, 0.0)
                )
            );

            zetaVib_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "zetav_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zetav_" + species_[i], dimless, 0.0)
                )
            );

            zetaElec_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "zetael_" + species_[i],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zetael_" + species_[i], dimless, 0.0)
                )
            );
        }
        else
        {
            if
            (
                not YHeader.typeHeaderOk<volScalarField>(false)
             && not XHeader.typeHeaderOk<volScalarField>(false)
            )
            {
                FatalErrorIn
                (
                    "basic2MultiComponentMixture::basic2MultiComponentMixture"
                )   << "Mass-fractions or molar-fractions header missing in "
                    << "the 0 folder" << nl;
            }
            if (not TtHeader.typeHeaderOk<volScalarField>(false))
            {
                FatalErrorIn
                (
                    "basic2MultiComponentMixture::basic2MultiComponentMixture"
                )   << "Translational temperature header missing in the 0/ "
                    << "folder" << nl;
            }
            FatalError<< exit(FatalError);
        }
    }

    /*forAll(spmodeTv_, speciei) // TODO ABORTIVE WORK
    {
        forAll(spmodeTv_[speciei], vibMode)
        {
            spmodeTv_[speciei][vibMode] = spTv_[speciei];
        }
    }*/

    e_ = new volScalarField
    (
        IOobject
        (
            "e",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("e", dimEnergy/dimMass, 0.0),
        this->he2BoundaryTypes
        (
            Tt.boundaryField(),
            downgradeToSingleTemperature
        ),
        this->he2BoundaryBaseTypes(Tt.boundaryField())
    );

    // Do not enforce constraint of sum of mass fractions to equal 1 here
    // - not applicable to all models
    {
        IOdictionary hTCPropertiesDict
        (
            IOobject
            (
                "hTCProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        writenD_ =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "numberDensity",
                false
            );
        writepD_ =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "partialDensity",
                false
            );
        writeX_ =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "molarFraction",
                false
            );
        writepP_ =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "partialPressure",
                false
            );
        writezeta_ =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "degreesOfFreedom",
                false
            );
        writehev_  =
            hTCPropertiesDict.subDict("mixtureOutputs").lookupOrDefault<bool>
            (
                "vibrationalEnergy",
                false
            );
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basic2MultiComponentMixture::write()
{
    /*forAll(spmodeTv_, speciei) // TODO ABORTIVE WORK
    {
        if (spmodeTv_[speciei].size() > 1)
        {
            forAll(spmodeTv_[speciei], vibMode)
            {
                spmodeTv_[speciei][vibMode].write();
            }
        }
    }*/

    if (writenD_)
    {
        forAll(species(), speciei)
        {
            nD_[speciei].write();
        }
    }

    if (writepD_)
    {
        forAll(species(), speciei)
        {
            pD_[speciei].write();
        }
    }

    if (writeX_)
    {
        forAll(species(), speciei)
        {
            X_[speciei].write();
        }
    }

    if (writepP_)
    {
        forAll(species(), speciei)
        {
            pP_[speciei].write();
        }
    }

    if (writezeta_)
    {
        forAll(species(), speciei)
        {
            zetaRot_[speciei].write();
            zetaVib_[speciei].write();
            zetaElec_[speciei].write();
        }
    }

    if (writehev_)
    {
        forAll(species(), speciei)
        {
            hev_[speciei].write();
            heel_[speciei].write();
            hevel_[speciei].write();
            e_->write();
        }
    }
}


// ************************************************************************* //
