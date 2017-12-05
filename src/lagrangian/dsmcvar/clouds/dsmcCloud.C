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


void Foam::dsmcCloud::buildCellOccupancyFromScratch()
{   
    cellOccupancy_.clear();
    cellOccupancy_.setSize(mesh_.nCells());

    buildCellOccupancy();
}

void Foam::dsmcCloud::buildCollisionSelectionRemainderFromScratch()
{       
    collisionSelectionRemainder_.clear();
    collisionSelectionRemainder_.setSize(mesh_.nCells());
    
    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, cO)
    {
        collisionSelectionRemainder_[cO] = rndGen_.scalar01();
    }
}

void Foam::dsmcCloud::resetBoundaries()
{  
    boundaries_.setNewConfig();
    boundaryMeas_.reset();
}

void Foam::dsmcCloud::resetMeasurementTools()
{   
    trackingInfo_.reset();
    fields_.resetFields();
    
    cellMeas_.reset(); // NEW VINCENT
}

void Foam::dsmcCloud::removeElectrons()
{   
    forAll(cellOccupancy_, c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];
        
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            
            const dsmcParcel::constantProperties& constProp 
                                = constProps(p->typeId());
                                
            const label& charge = constProp.charge();
            
            scalar RWF = 1.0;
            
            if(axisymmetric_)
            {
                const point& cC = mesh_.cellCentres()[c];
                scalar radius = cC.y();
                
                RWF = 1.0 + maxRWF_*(radius/radialExtent_);
            }
            
            scalar mass = constProps(p->typeId()).mass();
            
            momentumMean_[c] += mass*RWF*p->U();
            rhoMMean_[c] += mass*RWF;
            
            if(charge == -1)
            {
                rhoNMeanElectron_[c] += 1.0*RWF;
                rhoMMeanElectron_[c] += mass*RWF;
                momentumMeanElectron_[c] += mass*RWF*p->U();
                linearKEMeanElectron_[c] += mass*RWF*(p->U() & p->U());
                
                //found an electron
                deleteParticle(*p);
            }
        }
    }
}

void Foam::dsmcCloud::addElectrons()
{    
    
    label electronIndex = 0;
    
    label electronTypeId = -1;
            
    //find electron typeId
    forAll(constProps_, cP)
    {
        const label& electronCharge = constProps_[cP].charge();
        
        if(electronCharge == -1)
        {
            electronTypeId = cP;
            break;
        }
    }
    
    forAll(cellOccupancy_, c)
    {                
        if(rhoMMeanElectron_[c] > VSMALL)
        {
            
            scalar V = mesh_.cellVolumes()[c];
                
            scalar rhoMMeanElectron = rhoMMeanElectron_[c]*nParticle_[c]/V;
            scalar rhoNMeanElectron = rhoNMeanElectron_[c]*nParticle_[c]/V;
            vector UElectron = momentumMeanElectron_[c] /(rhoMMeanElectron*V);
            scalar linearKEMeanElectron = (0.5*linearKEMeanElectron_[c]*nParticle_[c])/V;
            
            electronTemperature_[c] = 2.0/(3.0*physicoChemical::k.value()*rhoNMeanElectron)
                                    *(linearKEMeanElectron - 0.5*rhoMMeanElectron*(UElectron & UElectron));
                                    
//             electronVelocity_[c] = UElectron;
        }
    }
        
    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();
        
        const dsmcParcel::constantProperties& constProp 
                            = constProps(p.typeId());
                            
        label charge = constProp.charge();
        
        if(charge  == 1)
        {
            //found an ion, add an electron here
            
            //electron temperature will be zero if there have been no
            //electrons in the cell during the simulation
            
            
            label cellI = p.cell();
            vector position = p.position();
            label tetFaceI = p.tetFace();
            label tetPtI = p.tetPt();
            
            if(electronTemperature_[cellI] < VSMALL)
            {
                electronTemperature_[cellI] = 6000.0;
            }
            if(electronTemperature_[cellI] > 8.0e4)
            {
                electronTemperature_[cellI] = 30000.0;
            }
                

            vector electronVelocity = equipartitionLinearVelocity
                (
                    electronTemperature_[cellI],
                    constProps_[electronTypeId].mass()
                );
              
            if(rhoMMean_[cellI] > VSMALL)
            {
                cellVelocity_[cellI] = momentumMean_[cellI]/rhoMMean_[cellI];
            }
            
            labelList vibLevel;
                
            electronVelocity += cellVelocity_[cellI];

            scalar RWF = p.RWF();

            addNewParcel
            (
                position,
                electronVelocity,
                RWF,
                0.0,
                0,
                cellI,
                tetFaceI,
                tetPtI,
                electronTypeId,
                0,
                0,
                vibLevel
            );
            
            electronIndex++;
        }            
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


void Foam::dsmcCloud::calculatenParticleField()
{
    // Find the cell with minimum volume
    // This cell will serve as a reference cell to set the 
    // nParticle/timeStep ratio
    
    label refCell = 0;
    scalar minVolume = Foam::GREAT;
    
    const scalarField& volumeCells = mesh_.V();
    
    forAll(nParticle_, cellI)
    {
        if(volumeCells[cellI] < minVolume)
        {
            minVolume = volumeCells[cellI];
            refCell = cellI;
        }
    }
    
    const scalar nParticleRef = nParticle_[refCell];
        
    forAll(nParticle_, cellI)
    {
        nParticle_[cellI] = nParticleRef*volumeCells[cellI]/minVolume;
    }
}


/*void Foam::dsmcCloud::setVariableTimeStepField(const label refCell)
{
    const scalar nParticleTimeStepRatio = nParticle_[refCell]
        /variableTimeStep_[refCell];
        
    forAll(variableTimeStep_, cellI)
    {
        if(cellI != refCell)
        {
            variableTimeStep_[cellI] = nParticle_[cellI]
                /nParticleTimeStepRatio;
        }
    }
}*/


void Foam::dsmcCloud::calculateTimeStepField()
{
    const scalarField& volumeCells = mesh_.V();
    const scalar minVolume = gMin(volumeCells);

    scalar nParticleTimeStepRatio = 0;
    label refCell = 0;
    
    forAll(nParticle_, cellI)
    {
        if(mag(volumeCells[cellI] - minVolume) < 1e-9)
        {
            nParticleTimeStepRatio = nParticle_[cellI]
                /variableTimeStep_[cellI];
            refCell = cellI;    
            break;
        }
    }
    
    if (Pstream::parRun())
    {
        reduce(nParticleTimeStepRatio, maxOp<scalar>());
    }
    
    forAll(variableTimeStep_, cellI)
    {
        variableTimeStep_[cellI] = nParticle_[cellI]
            /nParticleTimeStepRatio;
    }
    
    refQuantitiesVariableTimeStep_[0] = volumeCells[refCell];
    refQuantitiesVariableTimeStep_[1] = nParticle_[refCell];
    refQuantitiesVariableTimeStep_[2] = variableTimeStep_[refCell];
}


void Foam::dsmcCloud::updateVariableTimeStep()
{
    const scalar refnParticle = refQuantitiesVariableTimeStep_[1];
    const scalar refDeltaT = refQuantitiesVariableTimeStep_[2];
        
    if(variableTimeStepMethod_) // NEW VINCENT
    {
        const scalarField& volumeCells = mesh_.V();
        
        const scalar refVolume = refQuantitiesVariableTimeStep_[0];
        

        forAll(nParticle_, cellI)
        {
            nParticle_[cellI] = refnParticle*volumeCells[cellI]/refVolume;
        }
        
        const scalar nParticleTimeStepRatio = refnParticle/refDeltaT;
        
        forAll(variableTimeStep_, cellI)
        {
            variableTimeStep_[cellI] = nParticle_[cellI]
                /nParticleTimeStepRatio;
        }
    }
    else
    {
        forAll(nParticle_, cellI)
        {
            nParticle_[cellI] = refnParticle;
            variableTimeStep_[cellI] = refDeltaT;
        }
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


void Foam::dsmcCloud::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification,
    const labelList vibLevel
)
{
    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification,
        vibLevel
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
    const dynamicFvMesh& mesh,
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
    controlDict_
    (
        IOobject
        (
            "controlDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    variableTimeStepMethod_(particleProperties_.lookupOrDefault("variableTimeStepMethod", false)),
    nParticle_
    (
        IOobject
        (
            "nEquivalentParticles",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    variableTimeStep_
    (
        IOobject
        (
            "variableTimeStep",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("variableTimeStep", dimTime, mesh_.time().deltaTValue())
    ),
    refQuantitiesVariableTimeStep_(),
    axisymmetric_(Switch(particleProperties_.lookup("axisymmetricSimulation"))),
    radialExtent_(0.0),
    maxRWF_(1.0),
    nTerminalOutputs_(readLabel(controlDict_.lookup("nTerminalOutputs"))),
    cellOccupancy_(/*mesh_.nCells()*/), // EDITED VINCENT
    rhoNMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    momentumMeanElectron_(mesh_.nCells(), vector::zero),
    momentumMean_(mesh_.nCells(), vector::zero),
    linearKEMeanElectron_(mesh_.nCells(), 0.0),
    electronTemperature_(mesh_.nCells(), 0.0),
    cellVelocity_(mesh_.nCells(), vector::zero),
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
    collisionSelectionRemainder_(/*mesh_.nCells(), 0*/),
    constProps_(),
    rndGen_(label(clock::getTime()) + 7183*Pstream::myProcNo()), // different seed every time simulation is started - needed for ensemble averaging!
    controllers_(t, mesh, *this),
    dynamicLoadBalancing_(t, mesh, *this),
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
    cellMeas_(mesh, *this, true)
{
    if(variableTimeStepMethod_) // NEW VINCENT
    {
        calculateTimeStepField();
    }
    else
    {
        refQuantitiesVariableTimeStep_[1] = nParticle_[0];
        refQuantitiesVariableTimeStep_[2] = variableTimeStep_[0];
    }
    
    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();
    
    if(axisymmetric_)
    {
        radialExtent_ = readScalar(particleProperties_.lookup("radialExtentOfDomain"));
        maxRWF_ = readScalar(particleProperties_.lookup("maxRadialWeightingFactor"));
    }

    reactions_.initialConfiguration();

    //buildCellOccupancy(); // DELETED VINCENT
    buildCellOccupancyFromScratch(); // NEW VINCENT
    buildCollisionSelectionRemainderFromScratch(); // NEW VINCENT

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



// running dsmcInitialise
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const dynamicFvMesh& mesh,
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
    controlDict_
    (
        IOobject
        (
            "controlDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    variableTimeStepMethod_(particleProperties_.lookupOrDefault("variableTimeStepMethod", false)),
    nParticle_
    (
        IOobject
        (
            "nEquivalentParticles",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("nEquivalentParticles", dimless, 1.0)
    ),
    variableTimeStep_
    (
        IOobject
        (
            "variableTimeStep",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("variableTimeStep", dimTime, mesh_.time().deltaTValue())
    ),
    refQuantitiesVariableTimeStep_(),
    axisymmetric_(Switch(particleProperties_.lookup("axisymmetricSimulation"))),
    radialExtent_(0.0),
    maxRWF_(1.0),
    nTerminalOutputs_(readLabel(controlDict_.lookup("nTerminalOutputs"))),
    cellOccupancy_(),
    rhoNMeanElectron_(),
    rhoMMeanElectron_(),
    rhoMMean_(),
    momentumMeanElectron_(),
    momentumMean_(),
    linearKEMeanElectron_(),
    electronTemperature_(),
    cellVelocity_(),
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
    dynamicLoadBalancing_(t, mesh, *this),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    boundaryMeas_(mesh, *this),
    cellMeas_(mesh, *this)
{
    // NEW VINCENT ************************************************************
    const bool nParticleFromFile = particleProperties_.lookupOrDefault<bool>("nEquivalentParticlesFromFile", false);
    
    if(nParticleFromFile and variableTimeStepMethod_)
    {
        nParticle_.regIOobject::read();
    }
    else
    {
        const scalar nParticleFromDict = readScalar(particleProperties_.lookup("nEquivalentParticles"));
        
        forAll(nParticle_, cellI)
        {
            nParticle_[cellI] = nParticleFromDict;
        }
    }
    
    if(variableTimeStepMethod_)
    {
        calculatenParticleField();
    }
    // END NEW VINCENT ********************************************************
    
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
    
    if(axisymmetric_)
    {
        radialExtent_ = readScalar(particleProperties_.lookup("radialExtentOfDomain"));
        maxRWF_ = readScalar(particleProperties_.lookup("maxRadialWeightingFactor"));
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


Foam::dsmcCloud::~dsmcCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::dsmcCloud::evolve()
{
    boundaries_.updateTimeInfo();//****
    fields_.updateTimeInfo();//****
    controllers_.updateTimeInfo();//****

    dsmcParcel::trackingData td(*this);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();//****
    boundaries_.controlBeforeMove();//****
    
    //Remove electrons
    //removeElectrons();

    // Move the particles ballistically with their current velocities
    //Cloud<dsmcParcel>::move(td, mesh_.time().deltaTValue());
    Cloud<dsmcParcel>::move(td, deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    if(axisymmetric_)
    {
        axisymmetricWeighting();
        buildCellOccupancy();
    }

    //Add electrons back after the move function
    //addElectrons();
    
    // Update cell occupancy
    buildCellOccupancy();

    controllers_.controlBeforeCollisions();//****
    boundaries_.controlBeforeCollisions();//****

//     Info << "collisions" << endl;

    // Calculate new velocities via stochastic collisions
    collisions();

    controllers_.controlAfterCollisions();//****
    boundaries_.controlAfterCollisions();//****

    reactions_.outputData();

    fields_.calculateFields();//****
    fields_.writeFields();//****

    controllers_.calculateProps();//****
    controllers_.outputResults();//****

    boundaries_.calculateProps();//****
    boundaries_.outputResults();//****

    trackingInfo_.clean(); //****
    boundaryMeas_.clean(); //****
    cellMeas_.clean();  
}

Foam::label Foam::dsmcCloud::nTerminalOutputs()
{
    return nTerminalOutputs_;
}

void Foam::dsmcCloud::info() const
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    //scalar nMol = nDsmcParticles*nParticle_;
    scalar nMol = 0.0;
    forAll(cellOccupancy_, c)
    {
       nMol += cellOccupancy_[c].size()*nParticle_[c];
    }
    
    
    scalar linearKineticEnergy = infoMeasurements()[1];
    reduce(linearKineticEnergy, sumOp<scalar>());
    
    scalar rotationalEnergy = infoMeasurements()[2];
    reduce(rotationalEnergy, sumOp<scalar>());
    
    scalar vibrationalEnergy = infoMeasurements()[3];
    reduce(vibrationalEnergy, sumOp<scalar>());
    
    scalar electronicEnergy = infoMeasurements()[4];
    reduce(electronicEnergy, sumOp<scalar>());

//     vector linearMomentum = linearMomentumOfSystem();
//     reduce(linearMomentum, sumOp<vector>());

//     scalar linearKineticEnergy = linearKineticEnergyOfSystem();
//     reduce(linearKineticEnergy, sumOp<scalar>());
// 
//     scalar rotationalEnergy = rotationalEnergyOfSystem();
//     reduce(rotationalEnergy, sumOp<scalar>());
//     
//     scalar vibrationalEnergy = vibrationalEnergyOfSystem();
//     reduce(vibrationalEnergy, sumOp<scalar>());
//     
//     scalar electronicEnergy = electronicEnergyOfSystem();
//     reduce(electronicEnergy, sumOp<scalar>());

    Info << "    Number of dsmc particles        = "
    << nDsmcParticles
    << endl;

    if (nDsmcParticles)
    {
//         Info<< "    Number of molecules             = "
//             << nMol << nl
//             << "    Mass in system                  = "
//             << returnReduce(massInSystem(), sumOp<scalar>()) << nl
//             << "    Average linear momentum         = "
//             << linearMomentum/nMol << nl
//             << "    |Total linear momentum|         = "
//             << mag(linearMomentum) << nl
       Info << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = "
            << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = "
            << vibrationalEnergy/nMol << nl
            << "    Average electronic energy       = "
            << electronicEnergy/nMol << nl
            << "    Total energy                    = "
            << (rotationalEnergy + linearKineticEnergy
                + vibrationalEnergy + electronicEnergy)
            << endl;
    }
}

void Foam::dsmcCloud::loadBalanceCheck()
{
    dynamicLoadBalancing_.update();
}

void Foam::dsmcCloud::loadBalance(const int noRefinement)
{
    dynamicLoadBalancing_.perform(noRefinement);
}


void Foam::dsmcCloud::autoMap(const mapPolyMesh& mapper)
{    
    dsmcParcel::trackingData td(*this);

    //Cloud<dsmcParcel>::autoMap(td, mapper, true);
    Cloud<dsmcParcel>::autoMap(td, mapper);
    
    updateVariableTimeStep(); // NEW VINCENT

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();
    resetBoundaries();
    resetMeasurementTools();
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

Foam::labelList Foam::dsmcCloud::equipartitionVibrationalEnergyLevel
(
    scalar temperature,
    scalar vibrationalDof,
    label typeId
)
{
    labelList vibLevel(vibrationalDof, 0);

    if (vibrationalDof < SMALL)
    {
        return vibLevel;
    }
    else
    {  
        forAll(vibLevel, i)
        {
            label j = -log(rndGen_.scalar01())*temperature/constProps(typeId).thetaV()[i];
            vibLevel[i] = j;
        }
    }

    return vibLevel;
}

Foam::label Foam::dsmcCloud::equipartitionElectronicLevel
(
    scalar temperature,
    List<label> degeneracyList_,
    List<scalar> electronicEnergyList_,
    label typeId
)

{
    
    scalar EMax = physicoChemical::k.value()*temperature;
    label jMax = constProps(typeId).numberOfElectronicLevels();

    label jDash = 0; // random integer between 0 and jMax-1.
    scalar EJ = 0.0; // maximum possible electronic energy level within list based on k*TElec.
    label gJ = 0; // maximum possible degeneracy level within list.
    label jSelect = 0; // selected intermediate integer electronic level (0 to jMax-1).
    scalar expMax = 0.0; // maximum denominator value in Liechty pdf (see below).
    scalar expSum = 0.0; // Summation term based on random electronic level. 
    scalar boltz = 0; // Boltzmann distribution of Eq. 3.1.1 of Liechty thesis.
    scalar func = 0.0; // distribution function Eq. 3.1.2 of Liechty thesis.

        
    if (jMax == 1) // Only the electron E- has 1 degeneracy level = 0 (for programming purposes) in constant/dsmcProperties
    {   
//         return EEle;
        return 0;
    }
    if(temperature < VSMALL)
    {
        return jDash;
    }
    else
    {           
        // Calculate summation term in denominator of Eq.3.1.1 from Liechty thesis.     
        label i = 0;
        do
        {
            expSum += degeneracyList_[i]*exp((-electronicEnergyList_[i]/EMax));  
            i += 1;
        } while (i < jMax);     

        // select maximum integer energy level based on boltz value.
        // Note that this depends on the temperature.   
        
        scalar boltzMax = 0.0;
        
        for (label ii = 0; ii < jMax; ii++)
        {
            //Eq. 3.1.1 of Liechty thesis.   
            boltz = degeneracyList_[ii]*exp((-electronicEnergyList_[ii]/EMax))/expSum;
            
            if (boltzMax < boltz)
            {
                boltzMax = boltz;
                jSelect = ii;
            }               
        }

        EJ = electronicEnergyList_[jSelect]; //Max. poss energy in list : list goes from [0] to [jMax-1]
        gJ = degeneracyList_[jSelect]; //Max. poss degeneracy in list : list goes from [0] to [jMax-1]
        expMax = gJ*exp((-EJ/EMax)); // Max. in denominator of Liechty pdf for initialisation/
                                     //wall bcs/freestream EEle etc..
        
        do // acceptance - rejection based on Eq. 3.1.2 of Liechty thesis.      
        {               
            jDash = rndGen_.integer(0,jMax-1);
            func = degeneracyList_[jDash]*exp((-electronicEnergyList_[jDash]/EMax))/expMax;
        } while( !(func > rndGen_.scalar01()));
    }

    return jDash;
}

Foam::scalar Foam::dsmcCloud::postCollisionRotationalEnergy
(
    scalar rotationalDof,
    scalar ChiB  
)
{
    scalar energyRatio = 0.0;
    
    if(rotationalDof == 2.0)
    {
        energyRatio = 1.0 - pow(rndGen_.scalar01(),(1.0/ChiB));
    }
    else
    {
        scalar ChiA = 0.5*rotationalDof;
        
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
    }
    
    return energyRatio; 
}

Foam::label Foam::dsmcCloud::postCollisionVibrationalEnergyLevel
(
    bool postReaction,
    label vibLevel,
    label iMax,
    scalar thetaV,
    scalar thetaD,
    scalar refTempZv,
    scalar omega,
    scalar Zref,
    scalar Ec
)
{   
    label iDash = vibLevel;
    
    if(postReaction)
    {
        // post-collision quantum number
//         label iDash = 0; 
        scalar func = 0.0;
        scalar EVib = 0.0;

        do // acceptance - rejection 
        {
            iDash = rndGen_.integer(0,iMax);
            EVib = iDash*physicoChemical::k.value()*thetaV;
            
            // - equation 5.61, Bird
            func = pow((1.0 - (EVib / Ec)),(1.5 - omega));

        } while( !(func > rndGen_.scalar01()) );
    }
    else
    {
        // - "quantised collision temperature" (equation 3, Bird 2010), denominator from Bird 5.42

        scalar TColl = (iMax*thetaV) / (3.5 - omega);
        
        scalar pow1 = pow((thetaD/TColl),0.33333) - 1.0;

        scalar pow2 = pow ((thetaD/refTempZv),0.33333) -1.0;

        // - vibrational collision number (equation 2, Bird 2010)
        scalar ZvP1 = pow((thetaD/TColl),omega); 
        
        scalar ZvP2 = pow(Zref*(pow((thetaD/refTempZv),(-1.0*omega))),(pow1/pow2));
        
        scalar Zv = ZvP1*ZvP2;
        
        //In order to obtain the relaxation rate corresponding to Zv with the collision 
        //energy-based procedure, the inelastic fraction should be set to about 1/(5Zv)
        //Bird 2008 RGD "A Comparison of Collison Energy-Based and Temperature-Based..."
        
        scalar inverseVibrationalCollisionNumber = 1.0/(5.0*Zv);
        
//          scalar inverseVibrationalCollisionNumber = 1.0;

        if(inverseVibrationalCollisionNumber > rndGen_.scalar01())
        {
            // post-collision quantum number
            scalar func = 0.0;
            scalar EVib = 0.0;

            do // acceptance - rejection 
            {
                iDash = rndGen_.integer(0,iMax);
                EVib = iDash*physicoChemical::k.value()*thetaV;
                
                // - equation 5.61, Bird
                func = pow((1.0 - (EVib / Ec)),(1.5 - omega));

            } while( !(func > rndGen_.scalar01()) );
        }
    }
    
    return iDash;
}

Foam::label Foam::dsmcCloud::postCollisionElectronicEnergyLevel
(
    scalar Ec,
    label jMax,
    scalar omega,
    List<scalar> EElist,
    List<label> gList
)
{   
    label nPossStates = 0;
    label ELevel = -1;
        
    if(jMax == 1)
    {
        nPossStates = gList[0];
    }
    else
    {
        forAll(EElist, n)
        {
            if(Ec > EElist[n])
            {
                nPossStates += gList[n];
            }
        }
    }
    
    label II = 0;
    
    do
    {
        label nState = ceil(rndGen_.scalar01()*(nPossStates));
        label nAvailableStates = 0;
        label nLevel = -1;
        
        forAll(EElist, n)
        {
            nAvailableStates += gList[n];
            
            if(nState <= nAvailableStates && nLevel < 0)
            {
                nLevel = n;
            }
        }

        if(Ec > EElist[nLevel])
        {
            scalar prob = pow(1.0 - (EElist[nLevel]/Ec), 1.5-omega);
            
            if(prob > rndGen_.scalar01())
            {
                II = 1;
                ELevel = nLevel;
            }
        }
        
    }while(II==0);
    
    return ELevel;
        
//     label maxLev = 0;
//     label jSelectA = 0;
//     label jSelectB = 0;
//     label jSelect = 0;
//     label jDash = 0.0;
//     label gJ = 0.0;
//     label ELevel = 0;
//     scalar g = 0.0;
//     scalar gMax = 0.0;
//     scalar EJ = 0.0;
//     scalar denomMax = 0.0;
//     scalar prob = 0.0;
//     
//     // Determine the maximum possible integer energy level immediately below EcP.
//             
//     for (label ii = 0; ii < jMax; ii++)
//     {       
//         if (EElist[ii] > Ec)
//         {
//             break;
//         }
//         
//         maxLev = ii;
//         jSelectA = ii;
//         
//         //Eq. 3.1.6 of Liechty thesis.    
//         g = gList[ii]*pow((Ec - EElist[ii]),(1.5 - omega));
//         
//         if ( ii == 0 || gMax < g )
//         {
//             gMax = g;
//             jSelectB = ii;
//         }
//     }
// 
//     jSelect = jSelectA;
//     
//     if (jSelectB < jSelectA) 
//     {
//         jSelect = jSelectB; 
//     }
// 
//             
//     EJ = EElist[jSelect]; //Max. poss energy in list : list goes from [0] to [jSelect]
//     gJ = gList[jSelect]; //Max. poss degeneracy in list : list goes from [0] to [jSelect]
//     denomMax = gJ*pow((Ec - EJ), (1.5 - omega)); // Max. denominator of Liechty pdf for post-collision pdf.
//         
//     do // acceptance - rejection based on Eq. 3.1.8 of Liechty thesis.     
//     {               
//         jDash = rndGen_.integer(0,maxLev);
//         prob = gList[jDash]*pow((Ec - EElist[jDash]), (1.5 - omega))/denomMax;
// 
//     } while( !(prob > rndGen_.scalar01()));          
//     
//     ELevel = jDash;//post-collision Electronic energy.
//     
//     return ELevel;
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

void Foam::dsmcCloud::axisymmetricWeighting()
{
//     reBuildCellOccupancy();
    
    forAll(cellOccupancy_, c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
                        
            point cC = mesh_.cellCentres()[c];
            scalar radius = cC.y();
//             scalar radius = sqrt(sqr(p->position().y()) + sqr(p->position().z()));
            
            scalar oldRadialWeight = p->RWF();
            
            scalar newRadialWeight = 1.0;

            newRadialWeight = 1.0 + maxRWF_*(radius/radialExtent_);
            
            p->RWF() = newRadialWeight;
            
            if(oldRadialWeight > newRadialWeight) 
            {
                //particle might be cloned
                
                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;
                
                while(prob > 1.0)
                {
                    //add a particle and reduce prob by 1.0
                    
                    label cellI = p->cell();
                    vector position = p->position();
                    label tetFaceI = p->tetFace();
                    label tetPtI = p->tetPt();
                    
                    vector U = p->U();
                    
                    U.z() *= -1.0;

                    addNewParcel
                    (
                        position,
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        cellI,
                        tetFaceI,
                        tetPtI,
                        p->typeId(),
                        p->newParcel(),
                        p->classification(),
                        p->vibLevel()
                    );
                    
                    prob -= 1.0;
                }
                
                if(prob > rndGen_.scalar01())
                {
                    label cellI = p->cell();
                    vector position = p->position();
                    label tetFaceI = p->tetFace();
                    label tetPtI = p->tetPt();
                    
                    vector U = p->U();
                    
                    U.z() *= -1.0;

                    addNewParcel
                    (
                        position,
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        cellI,
                        tetFaceI,
                        tetPtI,
                        p->typeId(),
                        p->newParcel(),
                        p->classification(),
                        p->vibLevel()
                    );
                }
            }
            
            if(newRadialWeight > oldRadialWeight)
            {           
                //particle might be deleted
                if( (oldRadialWeight/newRadialWeight) < rndGen_.scalar01() )
                {
                    deleteParticle(*p);
                } 
            } 
        }
    }
    
    //reBuildCellOccupancy();
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
