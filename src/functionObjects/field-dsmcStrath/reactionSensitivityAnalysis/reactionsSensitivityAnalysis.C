/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "reactionsSensitivityAnalysis.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class chemistryType>
void Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
createFileNames()
{
    if (writeToFile() && !prodFilePtr_.valid())
    {
        prodFilePtr_ = createFile("production");
        writeHeader(prodFilePtr_(), "production");
        writeFileHeader(prodFilePtr_());

        consFilePtr_ = createFile("consumption");
        writeHeader(consFilePtr_(), "consumption");
        writeFileHeader(consFilePtr_());

        prodIntFilePtr_ = createFile("productionInt");
        writeHeader(prodIntFilePtr_(), "productionInt");
        writeFileHeader(prodIntFilePtr_());

        consIntFilePtr_ = createFile("consumptionInt");
        writeHeader(consIntFilePtr_(), "consumptionInt");
        writeFileHeader(consIntFilePtr_());
    }
}


template<class chemistryType>
void Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
writeFileHeader
(
    OFstream& os
)
{
    writeCommented(os, "Reaction");

    forAll(speciesNames_, k)
    {
        os << tab << speciesNames_[k] << tab;
    }

    os << nl << endl;
}


template<class chemistryType>
void Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
calculateSpeciesRR
(
    const basicChemistryModel& basicChemistry
)
{
    tmp<DimensionedField<scalar, volMesh>> RRt
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RR",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& RR = RRt.ref();

    scalar dt = time_.deltaT().value();

    endTime_ += dt;

    forAll(production_, speciei)
    {
        forAll(production_[speciei], reactioni)
        {
            RR = basicChemistry.calculateRR(reactioni, speciei);

            if (RR[0] > 0.0)
            {
                production_[speciei][reactioni] = RR[0];
                productionInt_[speciei][reactioni] =+ dt*RR[0];
            }
            else if (RR[0] < 0.0)
            {
                consumption_[speciei][reactioni] = RR[0];
                consumptionInt_[speciei][reactioni] =+ dt*RR[0];
            }
            else
            {
                production_[speciei][reactioni] = 0.0;
                consumption_[speciei][reactioni] = 0.0;
            }
        }
    }
}


template<class chemistryType>
void Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
writeSpeciesRR()
{

    consFilePtr_() << "time : " << mesh_.time().value() << tab << nl;
    consFilePtr_() << "delta T : "<< mesh_.time().deltaT().value() << nl << nl;
    prodFilePtr_() << "time : " << mesh_.time().value() << tab << nl;
    prodFilePtr_() << "delta T : "<< mesh_.time().deltaT().value() << nl << nl;

    consIntFilePtr_() << "start time : " << startTime_ << tab
        << "end time :" <<  endTime_ << nl;

    prodIntFilePtr_() << "start time : " << startTime_ << tab
        << "end time :" <<  endTime_ << nl;

    for (label reactioni = 0; reactioni < nReactions_; ++reactioni)
    {
        consFilePtr_() << reactioni << tab;
        consIntFilePtr_() << reactioni << tab;
        prodFilePtr_() << reactioni << tab;
        prodIntFilePtr_() << reactioni << tab;

        forAll(speciesNames_, i)
        {
            prodFilePtr_() << production_[i][reactioni] << tab;
            consFilePtr_() << consumption_[i][reactioni] << tab;
            prodIntFilePtr_() << productionInt_[i][reactioni] << tab;
            consIntFilePtr_() << consumptionInt_[i][reactioni] << tab;
            consumptionInt_[i][reactioni] = 0.0;
            productionInt_[i][reactioni] = 0.0;
        }
        consFilePtr_() << nl;
        consIntFilePtr_() << nl;
        prodFilePtr_() << nl;
        prodIntFilePtr_() << nl;
    }
    consFilePtr_() << nl << nl;
    consIntFilePtr_() << nl << nl;
    prodFilePtr_() << nl << nl;
    prodIntFilePtr_() << nl << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class chemistryType>
Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
reactionsSensitivityAnalysis
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    production_(0),
    consumption_(0),
    productionInt_(0),
    consumptionInt_(0),
    startTime_(0),
    endTime_(0),
    speciesNames_(),
    nReactions_(0),
    prodFilePtr_(),
    consFilePtr_(),
    prodIntFilePtr_(),
    consIntFilePtr_()
{
    read(dict);

    if (mesh_.nCells() != 1)
    {
        FatalErrorInFunction
            << "Function object only applicable to single cell cases"
            << abort(FatalError);
    }

    if (foundObject<basicChemistryModel>("chemistryProperties"))
    {
        const chemistryType& chemistry = refCast<const chemistryType>
        (
            lookupObject<basicChemistryModel>("chemistryProperties")
        );

        speciesNames_.setSize
        (
            chemistry.thermo().composition().species().size()
        );

        forAll(speciesNames_, i)
        {
            speciesNames_[i] = chemistry.thermo().composition().species()[i];
        }

        nReactions_ = chemistry.nReaction();

        if (production_.size() == 0)
        {
            production_.setSize(speciesNames_.size());
            consumption_.setSize(production_.size());
            productionInt_.setSize(production_.size());
            consumptionInt_.setSize(production_.size());

            forAll(production_, i)
            {
                production_[i].setSize(nReactions_, 0.0);
                consumption_[i].setSize(nReactions_, 0.0);
                productionInt_[i].setSize(nReactions_, 0.0);
                consumptionInt_[i].setSize(nReactions_, 0.0);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << " No chemistry model found. "
            << " Objects available are : " << mesh_.names()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class chemistryType>
Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
~reactionsSensitivityAnalysis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class chemistryType>
bool Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    return true;
}


template<class chemistryType>
bool Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
execute()
{
    createFileNames();

    const basicChemistryModel& chemistry =
        lookupObject<basicChemistryModel>("chemistryProperties");
    calculateSpeciesRR(chemistry);

    return true;
}


template<class chemistryType>
bool Foam::functionObjects::reactionsSensitivityAnalysis<chemistryType>::
write()
{
    if (Pstream::master())
    {
        writeSpeciesRR();

        startTime_ = endTime_;
    }

    return true;
}


// ************************************************************************* //
