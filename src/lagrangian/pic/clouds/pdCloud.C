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

#include "pdCloud.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

namespace Foam
{
//     defineParticleTypeNameAndDebug(pdParcel, 0);
    defineTemplateTypeNameAndDebug(Cloud<pdParcel>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::pdCloud::buildConstProps()
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

        constProps_[i] = pdParcel::constantProperties(molDict);
    }

}



void Foam::pdCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(pdCloud, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}

Foam::label Foam::pdCloud::pdkFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
    label entry = -1;
    label size = candidatesInCell.size();

    if(size > 0)
    {
        // choose a random number between 0 and the size of the candidateList size
        label randomIndex = rndGen_.position<label>(0, size - 1);
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

void Foam::pdCloud::updateCandidateSubList
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


Foam::label Foam::pdCloud::pdkFromCandidateSubList
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
        label randomIndex = rndGen_.position<label>(0, subCellSize - 1);
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

void Foam::pdCloud::collisions()
{
    collisionPartnerSelectionModel_->collide();
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


void Foam::pdCloud::addNewParcel
(
    const vector& position,
    const vector& U,
    /***********************************************************************/
    const vector& A,
    const scalar EPot,
    /***********************************************************************/
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
    pdParcel* pPtr = new pdParcel
    (
        mesh_,
        position,
        U,
        /***********************************************************************/
        A,
        EPot,
        /***********************************************************************/
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

Foam::scalar Foam::pdCloud::energyRatio
(
    scalar ChiA,
    scalar ChiB
)
{
    scalar ChiAMinusOne = ChiA - 1;

    scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
    {
        return rndGen_.sample01<scalar>();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen_.sample01<scalar>();

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
    } while (P < rndGen_.sample01<scalar>());

    return energyRatio;
}

Foam::scalar Foam::pdCloud::PSIm
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
        return rndGen_.sample01<scalar>();
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
        rPSIm = rndGen_.sample01<scalar>();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < rndGen_.sample01<scalar>());

    return rPSIm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// for running pdFOAM
Foam::pdCloud::pdCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<pdParcel>(mesh, cloudName, false),
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
    /***********************************************************************/
    emFields_(mesh, *this),
    electronModel_
	(
        ElectronModel::New
        (
            particleProperties_,
            *this
        )
    ),
    /***********************************************************************/
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
    boundaryMeas_(mesh, *this, true),
    measureDict_(particleProperties_.subDict(cloudName + "MeasurementOptions")),
    stdMeasure_(true),
    stdFrequency_(1),
    stdCount_(0),
    infoFrequency_(1)
{
    if (readFields)
    {
        pdParcel::readFields(*this);
    }

    //- check to see if non-standard measurement options have been selected

    if(measureDict_.found("standardFields"))
    {
        stdMeasure_     = readBool(measureDict_.lookup("standardFields"));
        if(stdMeasure_ && measureDict_.found("standardFieldFrequency"))
        {
            stdFrequency_   = readInt(measureDict_.lookup("standardFieldFrequency"));
        }
    }

    if(measureDict_.found("infoOutputFrequency"))
    {
        infoFrequency_   = readInt(measureDict_.lookup("infoOutputFrequency"));
    }


    buildConstProps();
//     Info << "Initial configuration" << endl;

    reactions_.initialConfiguration();

    buildCellOccupancy();


    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.sample01<scalar>();
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



// running pdIntialise
Foam::pdCloud::pdCloud
(
    Time& t,
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& pdInitialiseDict,
    const bool& clearFields
)
    :
    Cloud<pdParcel>(mesh, cloudName, false),
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
    /***********************************************************************/
    emFields_(mesh, *this, true),
    electronModel_(),
    /***********************************************************************/
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    boundaryMeas_(mesh, *this),
    measureDict_(particleProperties_.subDict(cloudName + "MeasurementOptions")),
    stdMeasure_(true),
    stdFrequency_(1),
    stdCount_(0),
    infoFrequency_(1)
{
    if(!clearFields)
    {
        pdParcel::readFields(*this);
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
    pdAllConfigurations conf(pdInitialiseDict, *this);
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


//Foam::pdCloud::~pdCloud()
//{Info << "pdCloud Destructor" << endl;}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pdCloud::evolve()
{
    /**************************************************************/
    /**              PIC/DSMC - Localisation and BCs             **/
    //cycles through models described in boundaryDict, fieldPropertiesDict and controllersDict
    //and updates their internal time stamp
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    stdCount_++;
    //- check to see if measuring using stdFields
    if(stdMeasure_ && stdCount_ == stdFrequency_)
    {
        // Reset the data collection for standard fields.
        standardFields_.resetFields();
    }

    emFields_.resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    /**************************************************************/
    /**               PD - Leapfrog method                      **/
    Info << "2. Electromagnetics" << endl;
    //- calculate fields based on particles at (x,v-0.5)_{k}
    //Info << "       - calculating fields" << endl;
    emFields_.calculateFields();

    //- calculate forces on part at x_{k}
    emFields_.calculateForces();
    //Info << "       - calculating lorentz force" << endl;

    Info << "3. Particle Integration" << endl;

    //Info << "       - leapfrog velocity" << endl;
    //- adjust velocity: (x,v-0.5) -> (x,v+0.5)
    pdParcel::trackingData td0(*this, 0);
    Cloud<pdParcel>::move(td0, mesh_.time().deltaTValue());

    //Info << "       - particle push" << endl;
    //- adjust position (ballistic part): (x,v+0.5) -> (x,v-0.5)
    pdParcel::trackingData td1(*this, 1);
    Cloud<pdParcel>::move(td1, mesh_.time().deltaTValue());

    //- build occupancy of particles
    buildCellOccupancy();

    /**************************************************************/
    /**                 DSMC - Collisions, Reactions             **/
    // Check for collisions and reactions at new location

    Info << "4. Collisions" << endl;

    // Apply controllers and boundary conditions before collisions
    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Update cell occupancy for reactions during collisions
    buildCellOccupancy();

    // Apply controller and boundary conditions after collisions
    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();

    /**************************************************************/
    /**        PIC/DSMC - Localisation, BC, Output, Cleanup      **/

    Info << "5. Writing" << endl;

    //- check to see if measuring using stdFields
    if(stdMeasure_ && stdCount_ == stdFrequency_)
    {
        // Calculate the volume field data
        standardFields_.calculateFields();
        stdCount_ = 0;
    }

    reactions_.outputData();

    fields_.calculateFields();
    fields_.writeFields();

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
}

void Foam::pdCloud::setupPIC()
{
    //cycles through models described in boundaryDict, fieldPropertiesDict and controllersDict
    //and updates their internal time stamp
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    //- reset the data collection for standard fields.
    standardFields_.resetFields();
    emFields_.resetFields();

    //- calculate initial fields and forces
    emFields().calculateFields();
    emFields().calculateForces();

    //- perform velocity adjust if required
    emFields().setupLeapfrog();
}


void Foam::pdCloud::info() const
{
    /** Initialise system variables **/
    const scalar e = electromagnetic::e.value();

    scalar nSpecies = typeIdList().size();

    List<vector> linearMomentum;
    List<vector> driftVelocity;

    List<scalar> potentialEnergy;
    List<scalar> rotationalEnergy;
    List<scalar> linearKineticEnergy;
    List<scalar> vibrationalEnergy;
    List<scalar> sysMass;
    List<scalar> sysCharge;

    List<scalar> nPart;

    linearMomentum.setSize(nSpecies,vector::zero);
    driftVelocity.setSize(nSpecies,vector::zero);

    potentialEnergy.setSize(nSpecies,0.0);
    rotationalEnergy.setSize(nSpecies,0.0);
    linearKineticEnergy.setSize(nSpecies,0.0);
    vibrationalEnergy.setSize(nSpecies,0.0);
    sysMass.setSize(nSpecies,0.0);
    sysCharge.setSize(nSpecies,0.0);

    nPart.setSize(nSpecies,0.0);

    // Number of pd parcels in system
    label nPdParticles = this->size();
    reduce(nPdParticles, sumOp<label>());

    /** Calculate species averages **/
    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();
        const scalar typeID = p.typeId();

        driftVelocity[typeID]  += p.U();
        nPart[typeID]          += 1.0;
    }

    //- Take averages
    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();
        const scalar typeID = p.typeId();

        driftVelocity[typeID] /= nPart[typeID];
    }

    //- Sum over proccessors
    reduce(nPart, sumOp<List<scalar> >());
    reduce(driftVelocity, sumOp<List<vector> >());

    /** Calculate system properties **/
    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();
        const scalar typeID = p.typeId();

        const pdParcel::constantProperties& cP = constProps
        (
            p.typeId()
        );

        linearMomentum[typeID]         += cP.mass()*p.U();

        linearKineticEnergy[typeID]    += 0.5*cP.mass()*(p.U() & p.U());
        potentialEnergy[typeID]        += 0.5*cP.Ze()*e*p.EPot();

        rotationalEnergy[typeID]       += p.ERot();
        vibrationalEnergy[typeID]      += p.EVib();

        sysMass[typeID]                += cP.mass();
        sysCharge[typeID]              += cP.Ze()*e;
    }

    /** Correct weighting **/
    forAll(nPart,tI)
    {
            nPart[tI]                  *= nParticle_;

            linearMomentum[tI]         *= nParticle_;
            linearKineticEnergy[tI]    *= nParticle_;
            potentialEnergy[tI]        *= nParticle_;
            rotationalEnergy[tI]       *= nParticle_;
            vibrationalEnergy[tI]      *= nParticle_;

            sysMass[tI]                *= nParticle_;
            sysCharge[tI]              *= nParticle_;
    }

    // Sum over processors
    reduce(linearMomentum, sumOp<List<vector> >());
    reduce(linearKineticEnergy, sumOp<List<scalar> >());

    reduce(potentialEnergy, sumOp<List<scalar> >());
    reduce(rotationalEnergy, sumOp<List<scalar> >());
    reduce(vibrationalEnergy, sumOp<List<scalar> >());

    reduce(sysCharge, sumOp<List<scalar> >());
    reduce(sysMass, sumOp<List<scalar> >());

    Info<< "Marker Cloud 4" << endl;

    /** Output data **/
    Info<< "Cloud name: " << this->name() << nl
        << "    Number of pd particles        = "
        << nPdParticles
        << endl;


    if (nPdParticles)
    {
        vector totalDriftVelocity           = vector::zero;
        vector totalLinearMomentum          = vector::zero;

        scalar totalSysMass                 = 0.0;
        scalar totalSysCharge               = 0.0;
        scalar totalPotentialEnergy         = 0.0;
        scalar totalLinearKineticEnergy     = 0.0;
        scalar totalRotationalEnergy        = 0.0;
        scalar totalVibrationalEnergy       = 0.0;

        // TODO: validate potential energy calculation for release version

        forAll(nPart,tI)
        {
            Info<< "Species: "
                << typeIdList()[tI] << nl
                << "    Number of " << typeIdList()[tI] << " molecules              = "
                << nPart[tI] << nl
                << "    Mass of " << typeIdList()[tI] << " in system                = "
                << sysMass[tI] << nl
                << "    Charge of " << typeIdList()[tI] << " system                 = "
                << sysCharge[tI] << nl
                //<< "    Average drift velocity of " << typeIdList()[tI] << "        = "
                //<< driftVelocity[tI] << nl
                //<< "    |Total linear momentum| of " << typeIdList()[tI] << "       = "
                //<< mag(linearMomentum[tI]) << nl
                //<< "    Potential energy of " << typeIdList()[tI] << "              = "
                //<< potentialEnergy[tI] << nl
                << "    Linear kinetic energy of " << typeIdList()[tI] << "   = "
                << linearKineticEnergy[tI] << nl
                << "    Rotational energy " << typeIdList()[tI] << "                = "
                << rotationalEnergy[tI] << nl
                << "    Vibrational energy " << typeIdList()[tI] << "               = "
                << vibrationalEnergy[tI] << nl
                << "    Total energy of " << typeIdList()[tI] << "                  = "
                //<< (rotationalEnergy[tI] + linearKineticEnergy[tI] + vibrationalEnergy[tI] + potentialEnergy[tI])
				<< (rotationalEnergy[tI] + linearKineticEnergy[tI] + vibrationalEnergy[tI])
                << endl;

                totalDriftVelocity           += driftVelocity[tI];
                totalLinearMomentum          += linearMomentum[tI];

                totalSysMass                 += sysMass[tI];
                totalSysCharge               += sysCharge[tI];
                totalPotentialEnergy         += potentialEnergy[tI];
                totalLinearKineticEnergy     += linearKineticEnergy[tI];
                totalRotationalEnergy        += rotationalEnergy[tI];
                totalVibrationalEnergy       += vibrationalEnergy[tI];
        }

        Info<< "Total: " << nl
                << "    Mass of system                          = "
                << totalSysMass << nl
                << "    Charge of system                        = "
                << totalSysCharge << nl
                //<< "    Average drift velocity of system        = "
                //<< totalDriftVelocity << nl
                << "    |Linear momentum| of system       = "
                << mag(totalLinearMomentum) << nl
                //<< "    Potential energy of particle system     = "
                //<< totalPotentialEnergy << nl
                << "    Linear kinetic energy of system   = "
                << totalLinearKineticEnergy << nl
                << "    Rotational energy of system             = "
                << totalRotationalEnergy << nl
                << "    Vibrational energy of system            = "
                << totalVibrationalEnergy << nl
                << "    Total energy of system                  = "
                //<< (totalRotationalEnergy + totalLinearKineticEnergy + totalVibrationalEnergy + totalPotentialEnergy)
				<< (totalRotationalEnergy + totalLinearKineticEnergy + totalVibrationalEnergy)
                << endl;
    }
}

Foam::vector Foam::pdCloud::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.GaussNormal<vector>();
}

Foam::scalar Foam::pdCloud::equipartitionRotationalEnergy
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
        ERot = -log(rndGen_.sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.sample01<scalar>());

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}

Foam::scalar Foam::pdCloud::equipartitionVibrationalEnergy
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
        label i = -log(rndGen_.sample01<scalar>())*temperature/constProps(typeId).thetaV();
        EVib = i*physicoChemical::k.value()*constProps(typeId).thetaV();
    }

    return EVib;
}

void Foam::pdCloud::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(pdCloud, *this, iter)
    {
        const pdParcel& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}

void Foam::pdCloud::reBuildCellOccupancy()
{
    buildCellOccupancy();
}

void Foam::pdCloud::insertParcelInCellOccupancy(pdParcel* p)
{
    cellOccupancy_[p->cell()].append(p);
    cellOccupancy_[p->cell()].shrink();
}

void Foam::pdCloud::removeParcelFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<pdParcel*> molsInCell(0);

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
