/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "dsmcCloud.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

namespace Foam
{
//     defineParticleTypeNameAndDebug(dsmcParcel, 0);
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::dsmcCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = dsmcParcel::constantProperties(molDict);
    }
}



void Foam::dsmcCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(dsmcCloud, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}

Foam::label Foam::dsmcCloud::pickFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
    label entry = -1;
    label size = candidatesInCell.size();

    if(size > 0)
    {
        // choose a random number between 0 and the size of the candidateList size
        label randomIndex = rndGen_.integer(0, size - 1);
        entry = candidatesInCell[randomIndex];

        // build a new list without the chosen entry
        DynamicList<label> newCandidates(0);
    
        forAll(candidatesInCell, i)
        {
            if(i != randomIndex)
            {
                newCandidates.append(candidatesInCell[i]);
            }
        }

        // transfer the new list
        candidatesInCell.transfer(newCandidates);
        candidatesInCell.shrink();
    }

    return entry;
}

void Foam::dsmcCloud::updateCandidateSubList
(
    const label& candidate,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " updating sub list (before) " << candidatesInSubCell << endl;

    label newIndex = findIndex(candidatesInSubCell, candidate);

    DynamicList<label> newCandidates(0);

    forAll(candidatesInSubCell, i)
    {
        if(i != newIndex)
        {
            newCandidates.append(candidatesInSubCell[i]);
        }
    }

    // transfer the new list
    candidatesInSubCell.transfer(newCandidates);
    candidatesInSubCell.shrink();

//     Info <<  " list (after) " << candidatesInSubCell << endl;
}


Foam::label Foam::dsmcCloud::pickFromCandidateSubList
(
    DynamicList<label>& candidatesInCell,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " list (before) " << candidatesInCell << endl;
//     Info << " sub list (before) " << candidatesInSubCell << endl;


    label entry = -1;
    label subCellSize = candidatesInSubCell.size();
    
    if(subCellSize > 0)
    {
        label randomIndex = rndGen_.integer(0, subCellSize - 1);
        entry = candidatesInSubCell[randomIndex];

//         Info<< "random index: " << randomIndex <<" entry " 
//             << entry << endl;

        DynamicList<label> newSubCellList(0);

        forAll(candidatesInSubCell, i)
        {
            if(i != randomIndex)
            {
                newSubCellList.append(candidatesInSubCell[i]);
            }
        }

        candidatesInSubCell.transfer(newSubCellList);
        candidatesInSubCell.shrink();

//         Info <<  " sub list (after) " << candidatesInSubCell << endl;

        label newIndex = findIndex(candidatesInCell, entry);

        DynamicList<label> newList(0);
    
        forAll(candidatesInCell, i)
        {
            if(i != newIndex)
            {
                newList.append(candidatesInCell[i]);
            }
        }

        candidatesInCell.transfer(newList);
        candidatesInCell.shrink();

//         Info <<  " list (after) " << candidatesInCell << endl;
    }

    return entry;
}

void Foam::dsmcCloud::collisions()
{
    collisionPartnerSelectionModel_->collide();
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


void Foam::dsmcCloud::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar ERot,
    const scalar EVib,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification
)
{
    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        ERot,
        EVib,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification
    );

    addParticle(pPtr);
}

Foam::scalar Foam::dsmcCloud::energyRatio
(
    scalar ChiA,
    scalar ChiB
)
{
    scalar ChiAMinusOne = ChiA - 1;

    scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
    {
        return rndGen_.scalar01();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen_.scalar01();

        if (ChiAMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiBMinusOne);
        }
        else if (ChiBMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiAMinusOne);
        }
        else
        {
            P =
                pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                    ChiAMinusOne
                )
               *pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                    /ChiBMinusOne,
                    ChiBMinusOne
                );
        }
    } while (P < rndGen_.scalar01());

    return energyRatio;
}

Foam::scalar Foam::dsmcCloud::PSIm
(
    scalar DOFm,
    scalar DOFtot
)
{
    if (DOFm == DOFtot)
    {
        return 1.0;
    }

    if (DOFm == 2.0 && DOFtot == 4.0)
    {
        return rndGen_.scalar01();
    }

    if (DOFtot < 4.0)
    {
        return (DOFm/DOFtot);
    }
    
    scalar rPSIm = 0.0;
    scalar prob = 0.0;

    scalar h1 = 0.5*DOFtot - 2.0;
    scalar h2 = 0.5*DOFm - 1.0 + 1.0e-5;
    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + 1.0e-5;

    do
    {
        rPSIm = rndGen_.scalar01();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < rndGen_.scalar01());

    return rPSIm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// for running dsmcFOAM
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(mesh_.nCells(), 0),
    constProps_(),
//     rndGen_(label(149382906) + 7183*Pstream::myProcNo()),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo()), // different seed every time simulation is started - needed for ensemble averaging!
    controllers_(t, mesh, *this),
    standardFields_(mesh, *this),
    fields_(t, mesh, *this),
    boundaries_(t, mesh, *this),
    trackingInfo_(mesh, *this, true),
    binaryCollisionModel_
    (
        BinaryCollisionModel::New
        (
            particleProperties_,
            *this
        )
    ),
	collisionPartnerSelectionModel_(),
    reactions_(t, mesh, *this),
    boundaryMeas_(mesh, *this, true)
{
    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();
//     Info << "Initial configuration" << endl;

    reactions_.initialConfiguration();

    buildCellOccupancy();

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.scalar01();
    }
    
	collisionPartnerSelectionModel_ = autoPtr<collisionPartnerSelection>
	(
    	collisionPartnerSelection::New(mesh, *this, particleProperties_)
	);

	collisionPartnerSelectionModel_->initialConfiguration();

    fields_.createFields();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
}



// running dsmcIntialise
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict,
    const bool& clearFields
)
    :
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 3, -1, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_(),
    constProps_(),
//     rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    rndGen_(label(clock::getTime()) + 1526*Pstream::myProcNo()), // different seed every time simulation is started - needed for ensemble averaging!
    controllers_(t, mesh),
    standardFields_(mesh, *this, true),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    boundaryMeas_(mesh, *this)
{
    if(!clearFields)
    {
        dsmcParcel::readFields(*this);
    }

    label initialParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(initialParcels, sumOp<label>());
    }

    if(clearFields)
    {
        Info << "clearing existing field of parcels " << endl;

        clear();

        initialParcels = 0;
    }

    buildConstProps();
    dsmcAllConfigurations conf(dsmcInitialiseDict, *this);
    conf.setInitialConfig();

    label finalParcels = this->size();
    
    if (Pstream::parRun())
    {
        reduce(finalParcels, sumOp<label>());
    }

    Info << nl << "Initial no. of parcels: " << initialParcels 
         << " added parcels: " << finalParcels - initialParcels
         << ", total no. of parcels: " << finalParcels 
         << endl;

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// Foam::dsmcCloud::~dsmcCloud()
// {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dsmcCloud::evolve()
{
    boundaries_.updateTimeInfo();//****
    fields_.updateTimeInfo();//****
    controllers_.updateTimeInfo();//****

    dsmcParcel::trackingData td(*this);

    // Reset the data collection fields
    standardFields_.resetFields();
//     resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();//****
    boundaries_.controlBeforeMove();//****

    // Move the particles ballistically with their current velocities
    Cloud<dsmcParcel>::move(td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    controllers_.controlBeforeCollisions();//****
    boundaries_.controlBeforeCollisions();//****
    Info << "collisions" << endl;

    // Calculate new velocities via stochastic collisions
    collisions();

    buildCellOccupancy(); //*** (for reactions)

    controllers_.controlAfterCollisions();//****
    boundaries_.controlAfterCollisions();//****

    // Calculate the volume field data
    standardFields_.calculateFields();

    reactions_.outputData();

    fields_.calculateFields();//****
//    fields_.writeFields();//****

    controllers_.calculateProps();//****
//    controllers_.outputResults();//****

    boundaries_.calculateProps();//****
//    boundaries_.outputResults();//****

    trackingInfo_.clean(); //****
    boundaryMeas_.clean(); //****
}

void Foam::dsmcCloud::evolveWrite()
{
    boundaries_.updateTimeInfo();//****
    fields_.updateTimeInfo();//****
    controllers_.updateTimeInfo();//****

    dsmcParcel::trackingData td(*this);

    // Reset the data collection fields
    standardFields_.resetFields();
//     resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();//****
    boundaries_.controlBeforeMove();//****

    // Move the particles ballistically with their current velocities
    Cloud<dsmcParcel>::move(td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    controllers_.controlBeforeCollisions();//****
    boundaries_.controlBeforeCollisions();//****
    Info << "collisions" << endl;

    // Calculate new velocities via stochastic collisions
    collisions();

    buildCellOccupancy(); //*** (for reactions)

    controllers_.controlAfterCollisions();//****
    boundaries_.controlAfterCollisions();//****

    // Calculate the volume field data
    standardFields_.calculateFields();

    reactions_.outputData();

    fields_.calculateFields();//****
//    fields_.writeFields();//****

    controllers_.calculateProps();//****
//    controllers_.outputResults();//****

    boundaries_.calculateProps();//****
//    boundaries_.outputResults();//****

    trackingInfo_.clean(); //****
    boundaryMeas_.clean(); //****
    
    List<label> cellsToWrite_;
    cellsToWrite_.setSize(3);
    cellsToWrite_[0] = 2;
    cellsToWrite_[1] = 247;
    cellsToWrite_[2] = 497;
    forAll(cellsToWrite_, cellI)
    {
        if(mesh_.nCells() > cellsToWrite_[cellI] - 1)
        {
            forAll(cellOccupancy_[cellsToWrite_[cellI]], pI)
            {
                dsmcParcel* p = cellOccupancy_[cellsToWrite_[cellI]][pI];
                //"identifier,u,v,w,eRot,eVib"
                Info << "miau" << cellI << "," << p->U()[1] << ","
                    << p->ERot() << "," << p->EVib() << endl;
            }
        }
    }
}




void Foam::dsmcCloud::info() const
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    scalar nMol = nDsmcParticles*nParticle_;

    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalEnergy = rotationalEnergyOfSystem();
    reduce(rotationalEnergy, sumOp<scalar>());
    
    scalar vibrationalEnergy = vibrationalEnergyOfSystem();
    reduce(vibrationalEnergy, sumOp<scalar>());

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = "
        << nDsmcParticles
        << endl;

    if (nDsmcParticles)
    {
        Info<< "    Number of molecules             = "
            << nMol << nl
            << "    Mass in system                  = "
            << returnReduce(massInSystem(), sumOp<scalar>()) << nl
            << "    Average linear momentum         = "
            << linearMomentum/nMol << nl
            << "    |Total linear momentum|         = "
            << mag(linearMomentum) << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = "
            << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = "
            << vibrationalEnergy/nMol << nl
            << "    Average total energy            = "
            << (rotationalEnergy + linearKineticEnergy + vibrationalEnergy)/*/nMol*/
            << endl;
    }
}


void Foam::dsmcCloud::resetHybrid
(
    volScalarField& TtrInitial,
    volVectorField& UInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qInitial,
    PtrList<volTensorField>& tauInitial,
    dimensionedScalar& B,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
    
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtrInitial[cellI];
                        scalar rotationalTemperature = TtrInitial[cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocityMiu
                            (
                                translationalTemperature,
                                cP.mass(),
                                B.value(),
                                qInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );
        
                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.vibrationalDegreesOfFreedom(),
                                i
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}

void Foam::dsmcCloud::resetHybrid2
(
    volScalarField& TtrInitial,
    volVectorField& UInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qInitial,
    PtrList<volTensorField>& tauInitial,
    dimensionedScalar& B,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
    
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtrInitial[cellI];
                        scalar rotationalTemperature = TtrInitial[cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocity
                            (
                                translationalTemperature,
                                cP.mass(),
                                qInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );
        
                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.vibrationalDegreesOfFreedom(),
                                i
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}

// Hybrid Reset----------------------------------------------------------------
void Foam::dsmcCloud::resetHybridMax
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& numberDensitiesField,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    for(label zoneToReset = 0; zoneToReset < zonesToReset.size(); zoneToReset++)
    {
        word regionName(zonesToReset[zoneToReset]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybridChapEnsk")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
    
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = 0;
                        scalar vibrationalTemperature = 0;
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->equipartitionLinearVelocity
                            (
                                translationalTemperature,
                                cP.mass()
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );
        
                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.vibrationalDegreesOfFreedom(),
                                i
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}


void Foam::dsmcCloud::resetHybridTra
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qtInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    for(label zoneToReset = 0; zoneToReset < zonesToReset.size(); zoneToReset++)
    {
        word regionName(zonesToReset[zoneToReset]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybridChapEnsk")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
    
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = 0;
                        scalar vibrationalTemperature = 0;
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocity
                            (
                                translationalTemperature,
                                cP.mass(),
                                qtInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );
        
                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.vibrationalDegreesOfFreedom(),
                                i
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}

void Foam::dsmcCloud::resetHybridTraRotVib
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    scalar time1_ = mesh_.time().elapsedCpuTime();
    molsToDeleteHybrid(mesh_, *this, typeOfReset);
    scalar time2_ = mesh_.time().elapsedCpuTime();
    Pout << "Proc " << UPstream::myProcNo() << " deletion time: " << time2_ - time1_ << "s (" << time1_ << ", " << time2_ << ")" << nl << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        List<scalar> massToIntroduce(TvInitial.size(), 0.0);//////////////
        List<scalar> massIntroduced(TvInitial.size(), 0.0);///////////////
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
                forAll(typeIdList_, i)  ////////////////////////////////
                {                       ////////////////////////////////
                    massToIntroduce[i] += numberDensitiesField[i][cellI] * mesh_.V()[cellI];
                }                       ////////////////////////////////
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }
                        massIntroduced[i] += nParticlesToInsert;///////////////
                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U;

                            scalar ERot;
                            scalar EVib;

                            this->generalisedChapmanEnskog
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                cP.mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                EVib,
                                U
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
        Info << "For zone " + regionName + ":" << endl;
        forAll(typeIdList_, i)
        {
            scalar mass = this->constProps(i).mass();
            massToIntroduce[i] *= mass * nParticle_;
            massIntroduced[i] *= mass * nParticle_;
            Info << "    For specie " << i << ", mTI: "
                << massToIntroduce[i] << "; mI: " << massIntroduced[i]
                << nl << "        (diff of " << 100.0 * (massIntroduced[i]
                / (massToIntroduce[i] + VSMALL) - 1.0) << "%)" << endl;
        }
    }
    buildCellOccupancy();
}

void Foam::dsmcCloud::resetHybridTraRotVib2
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if(zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];
    
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );
    
                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.scalar01()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U;

                            scalar ERot;
                            scalar EVib;

                            this->generalisedChapmanEnskog2
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                cP.mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                EVib,
                                U
                            );

                            U += velocity;
                
                            label newParcel = 0;
            
                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}
//-----------------------------------------------------------------------------

void Foam::dsmcCloud::shockReset()
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    const IOdictionary& shockDict
    (
        IOobject
        (
            "shockDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar initialParticles = readScalar(shockDict.lookup("initialParticles"));
    scalar maxPercentageParticles = readScalar(shockDict.lookup("maxPercentageParticles"));
    scalar rhoLeft = readScalar(shockDict.lookup("rhoLeft"));
    scalar rhoRight = readScalar(shockDict.lookup("rhoRight"));
    scalar xLeft = readScalar(shockDict.lookup("xLeft"));
    scalar xRight = readScalar(shockDict.lookup("xRight"));
    scalar areaTube = readScalar(shockDict.lookup("areaTube"));

    scalar maxDeltaParticles = maxPercentageParticles * initialParticles / 100.0;
 
    scalar deltaParticles = nDsmcParticles - initialParticles;
    Info << "deltaParticles [%]: " << 100.0 * deltaParticles / initialParticles
        << nl << "updating:" << (mag(deltaParticles) >= maxDeltaParticles)
        << endl;
    scalar deltaX = deltaParticles * nParticle_ * this->constProps(0).mass()
        / areaTube / (rhoRight - rhoLeft);
    Info << "deltaX: " << deltaX << endl;

    if(mag(deltaParticles) >= maxDeltaParticles)
    {

        forAll(mesh_.cells(), cellI)
        {

            const List<dsmcParcel*>& parcelsInCell
                = cellOccupancy_[cellI];

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                // deltaX > 0: if position > xRight + deltaX => out of domain
                // deltaX < 0: if position < xLeft + deltaX => out of domain
                if
                (
                    (p->classification() < 1000)
                    &&
                    (
                        (p->position().x() < xLeft + deltaX)
                            || (p->position().x() > xRight + deltaX)
                    )
                )
                {
                    p->classification() += 2000;
                    dsmcParcel* pPtr = new dsmcParcel(*p);
                    addParticle(pPtr);
                    p->classification() -= 2000;
                }
            }
        }
        forAllIter(dsmcCloud, *this, iter)
        {
            dsmcParcel& p = iter();
            if (p.classification() > 1000)
            {
                p.classification() -= 2000;
            }
            else
            {
                p.position().x() += deltaX;
                if((p.position().x() > xRight) || (p.position().x() < xLeft))
                {
                    deleteParticle(p);
                }
            }
        }
        buildCellOccupancy();
    }
}

Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *vector
        (
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal()
        );
}

Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    scalar temperature,
    scalar rotationalDof
)
{
    scalar ERot = 0.0;

    if (rotationalDof < SMALL)
    {
        return ERot;
    }
    else if (rotationalDof < 2.0 + SMALL && rotationalDof > 2.0 - SMALL)
    {
        // Special case for rDof = 2, i.e. diatomics;
        ERot = -log(rndGen_.scalar01())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.scalar01();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.scalar01());

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}

Foam::scalar Foam::dsmcCloud::equipartitionVibrationalEnergy
(
    scalar temperature,
    scalar vibrationalDof,
    label typeId
)
{
    scalar EVib = 0.0;

    if (vibrationalDof < SMALL)
    {
        return EVib;
    }
    else
    {  
        label i = -log(rndGen_.scalar01())*temperature/constProps(typeId).thetaV();
        EVib = i*physicoChemical::k.value()*constProps(typeId).thetaV();
    }

    return EVib;
}

Foam::scalar Foam::dsmcCloud::equipartitionVibrationalEnergy2
(
    scalar temperature,
    scalar vibrationalDof,
    label typeId
)
{
    scalar EVib = 0.0;

    if (vibrationalDof < SMALL)
    {
        return EVib;
    }
    else
    {  
        label i = -log(rndGen_.scalar01())*temperature/constProps(typeId).thetaV();
        Info << "Real: " <<
            -log(rndGen_.scalar01())*temperature/constProps(typeId).thetaV()
            << ", i: " << i << endl;
        EVib = i*physicoChemical::k.value()*constProps(typeId).thetaV();
    }

    return EVib;
}


Foam::vector Foam::dsmcCloud::chapmanEnskogVelocity
(
    scalar temperature,
    scalar mass,
    vector q,
    tensor tau
)
{
    scalar B = max(mag(q), mag(tau));
    scalar A = 1.0 + 30.0 * B;

    bool repeatTry = true;
    vector CTry;
    while (repeatTry)
    {
        CTry = vector
            (
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal()
            ) / sqrt(2.0);
        scalar gammaTry = 1.0 + (q & CTry) * (0.4 * (CTry & CTry) - 1.0)
            - (CTry & (tau & CTry));
        scalar rn = rndGen_.scalar01();
        if (gammaTry >= A * rn) repeatTry = false;
    }
    return CTry * sqrt(2.0 * physicoChemical::k.value() * temperature / mass);
}


Foam::vector Foam::dsmcCloud::chapmanEnskogVelocityMiu
(
    scalar temperature,
    scalar mass,
    scalar B,
    vector q,
    tensor tau
)
{
    scalar A = 1.0 + 30.0 * B;
    bool repeatTry = true;
    vector CTry;
    while (repeatTry)
    {
        CTry = vector
            (
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal()
            ) / sqrt(2.0);
        scalar gammaTry = 1.0 + (q & CTry) * (0.4 * (CTry & CTry) - 1.0)
            - (CTry & (tau & CTry));
        scalar rn = rndGen_.scalar01();
        if (gammaTry >= A * rn) repeatTry = false;
    }
    return CTry * sqrt(2.0 * physicoChemical::k.value() * temperature / mass);
}

void Foam::dsmcCloud::generalisedChapmanEnskog
(
    label typeID,
    scalar translationalTemperature,
    scalar rotationalTemperature,
    scalar vibrationalTemperature,
    scalar mass,
    vector D,
    vector qTra,
    vector qRot,
    vector qVib,
    tensor tau,
    scalar& ERot,
    scalar& EVib,
    vector& U
)
{
    scalar k = physicoChemical::k.value();
    scalar vibDOF = constProps(typeID).vibrationalDegreesOfFreedom();

    scalar B = max(mag(D), mag(tau));
    B = max(B, mag(qTra));
    B = max(B, mag(qRot));
    B = max(B, mag(qVib));
    scalar A = 1.0 + 30.0 * B;

    scalar thetaV = constProps(typeID).thetaV();
    scalar epsRotAv = rotationalTemperature / translationalTemperature;
    scalar epsVibAv = 0.0;
    if(vibDOF > 0 && vibrationalTemperature > SMALL)
    {
        epsVibAv = thetaV / vibrationalTemperature / (exp(thetaV
        / vibrationalTemperature) - 1.0);
    }

    vector CTry;
    scalar epsRot = 0.0;
    scalar epsVib = 0.0;
    bool repeatTry = true;
    while (repeatTry)
    {
        ERot = this->equipartitionRotationalEnergy
            (
                rotationalTemperature,
                constProps(typeID).rotationalDegreesOfFreedom()
            );
        
        EVib = this->equipartitionVibrationalEnergy
            (
                vibrationalTemperature,
                vibDOF,
                typeID
            );


        epsRot = ERot / k / translationalTemperature;
        if(vibDOF > 0 && vibrationalTemperature > SMALL)
        {
            epsVib = EVib / k / vibrationalTemperature;
        }

        CTry = vector
            (
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal()
            ) / sqrt(2.0);
        scalar gammaTry = 1.0 + 2.0 * (D & CTry) + (qTra & CTry) * (0.4
            * (CTry & CTry) - 1.0) + (qRot & CTry) * (epsRot - epsRotAv)
            + (qVib & CTry) * (epsVib - epsVibAv) - (CTry & (tau & CTry));
        scalar rn = rndGen_.scalar01();
        if (gammaTry >= A * rn) repeatTry = false;
    }

    U = CTry * sqrt(2.0 * k * translationalTemperature / mass);
}

void Foam::dsmcCloud::generalisedChapmanEnskog2
(
    label typeID,
    scalar translationalTemperature,
    scalar rotationalTemperature,
    scalar vibrationalTemperature,
    scalar mass,
    vector D,
    vector qTra,
    vector qRot,
    vector qVib,
    tensor tau,
    scalar& ERot,
    scalar& EVib,
    vector& U
)
{
    scalar B = max(mag(D), mag(tau));
    B = max(B, mag(qTra));
    B = max(B, mag(qRot));
    B = max(B, mag(qVib));
    scalar A = 1.0 + 30.0 * B;

    scalar thetaV = constProps(typeID).thetaV();
    scalar epsRotAv = rotationalTemperature / translationalTemperature;
    scalar epsVibAv = 0.0;
    if(vibrationalTemperature > SMALL)
    {
        epsVibAv = thetaV / vibrationalTemperature / (exp(thetaV
        / vibrationalTemperature) - 1.0);
    }

    vector CTry;
    scalar gammaTry = 1.0;
    scalar rn = 1.0;
    scalar epsRot = 0.0;
    scalar epsVib = 0.0;
    ERot = this->equipartitionRotationalEnergy
        (
            rotationalTemperature,
            constProps(typeID).rotationalDegreesOfFreedom()
        );
        
    EVib = this->equipartitionVibrationalEnergy
        (
            vibrationalTemperature,
            constProps(typeID).vibrationalDegreesOfFreedom(),
            typeID
        );

    epsRot = ERot / physicoChemical::k.value()
        / translationalTemperature;
    if(vibrationalTemperature > SMALL)
    {
        epsVib = EVib / physicoChemical::k.value()
            / vibrationalTemperature;
    }
    do
    {
        CTry = vector
            (
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal()
            ) / sqrt(2.0);
        gammaTry = 1.0 + 2.0 * (D & CTry) + (qTra & CTry) * (0.4
            * (CTry & CTry) - 1.0) + (qRot & CTry) * (epsRot - epsRotAv)
            + (qVib & CTry) * (epsVib - epsVibAv) - (CTry & (tau & CTry));
        rn = rndGen_.scalar01();
    } while (gammaTry >= A * rn);

    U = CTry * sqrt(2.0 * physicoChemical::k.value() * translationalTemperature
        / mass);
}

void Foam::dsmcCloud::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}

void Foam::dsmcCloud::reBuildCellOccupancy()
{
    buildCellOccupancy();
}

void Foam::dsmcCloud::insertParcelInCellOccupancy(dsmcParcel* p)
{
    cellOccupancy_[p->cell()].append(p);
    cellOccupancy_[p->cell()].shrink();
}

void Foam::dsmcCloud::removeParcelFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<dsmcParcel*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if(c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    molsInCell.shrink();
    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}


// ************************************************************************* //
