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
    boundaryMeas_.reset();
    boundaries_.setNewConfig();
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
            
            const dsmcParcel::constantProperties& constProp =
                constProps(p->typeId());
                                
            const label& charge = constProp.charge();
            
            const scalar& RWF = getRWF_cell(c);
            
            const scalar mass = constProps(p->typeId()).mass();
            
            momentumMean_[c] += mass*RWF*p->U();
            rhoMMean_[c] += mass*RWF;
            
            if(charge == -1)
            {
                rhoNMeanElectron_[c] += 1.0*RWF;
                rhoMMeanElectron_[c] += mass*RWF;
                momentumMeanElectron_[c] += mass*RWF*p->U();
                linearKEMeanElectron_[c] += mass*RWF*(p->U() & p->U());
                
                //- found an electron
                deleteParticle(*p);
            }
        }
    }
}

void Foam::dsmcCloud::addElectrons()
{      
    label electronTypeId = -1;
            
    //- find electron typeId
    forAll(constProps_, cP)
    {
        const label& particleCharge = constProps_[cP].charge();
        
        if(particleCharge == -1)
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
                
            scalar rhoMMeanElectron = rhoMMeanElectron_[c]*nParticle_/V;
            scalar rhoNMeanElectron = rhoNMeanElectron_[c]*nParticle_/V;
            vector UElectron = momentumMeanElectron_[c] /(rhoMMeanElectron*V);
            scalar linearKEMeanElectron = 
		            (0.5*linearKEMeanElectron_[c]*nParticle_)/V;
            
            electronTemperature_[c] = 2.0/(3.0*physicoChemical::k.value()
                * rhoNMeanElectron)*(linearKEMeanElectron 
                - 0.5*rhoMMeanElectron
                * (UElectron & UElectron));
        }
        
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];
        
        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            
            const dsmcParcel::constantProperties& constProp 
                                = constProps(p->typeId());
                                
            label charge = constProp.charge();
        
            if(charge == 1)
            {
                const label& cellI = p->cell();
                
                //- found an ion, add an electron here
                
                //- electron temperature will be zero if there have been no
                //  electrons in the cell during the simulation
                
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
                    cellVelocity_[cellI] = momentumMean_[cellI]
                        /rhoMMean_[cellI];
                }
                
                labelList vibLevel;
                    
                electronVelocity += cellVelocity_[cellI];

                addNewParcel
                (
                    p->position(),
                    electronVelocity,
                    p->RWF(),
                    0.0,
                    0,
                    p->cell(),
                    p->tetFace(),
                    p->tetPt(),
                    electronTypeId,
                    -1,
                    0,
                    vibLevel
                );
            }
        }
    }
        
//     forAllConstIter(dsmcCloud, *this, iter)
//     {
//         const dsmcParcel& p = iter();
//         
//         const dsmcParcel::constantProperties& constProp 
//                             = constProps(p.typeId());
//                             
//         label charge = constProp.charge();
//         
//         if(charge  == 1)
//         {
//             //found an ion, add an electron here
//             
//             //electron temperature will be zero if there have been no
//             //electrons in the cell during the simulation
//             
//             
//             label cellI = p.cell();
//             vector position = p.position();
//             label tetFaceI = p.tetFace();
//             label tetPtI = p.tetPt();
//             
//             if(electronTemperature_[cellI] < VSMALL)
//             {
//                 electronTemperature_[cellI] = 6000.0;
//             }
//             if(electronTemperature_[cellI] > 8.0e4)
//             {
//                 electronTemperature_[cellI] = 30000.0;
//             }
//                 
// 
//             vector electronVelocity = equipartitionLinearVelocity
//                 (
//                     electronTemperature_[cellI],
//                     constProps_[electronTypeId].mass()
//                 );
//               
//             if(rhoMMean_[cellI] > VSMALL)
//             {
//                 cellVelocity_[cellI] = momentumMean_[cellI]/rhoMMean_[cellI];
//             }
//             
//             labelList vibLevel;
//                 
//             electronVelocity += cellVelocity_[cellI];
// 
//             scalar RWF = p.RWF();
//             
//             addNewParcel
//             (
//                 position,
//                 electronVelocity,
//                 RWF,
//                 0.0,
//                 0,
//                 cellI,
//                 tetFaceI,
//                 tetPtI,
//                 electronTypeId,
//                 -1,
//                 0,
//                 vibLevel
//             );
//         }            
//     }
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
        // choose a random number between 0 and the size of the candidateList
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
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification,
    const labelList& vibLevel
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
    
    if(measureMediumProperties_ and newParcel != -1)
    {
        if(trackingProbability_ > rndGen_.scalar01())
        {
            pPtr->setTracked
            (
                true, 
                newParcel,
                mesh_.time().value(), 
                pPtr->position()
            );
        }
    }
    
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

// for running dsmcFoam+
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
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    axisymmetric_(particleProperties_.lookupOrDefault<Switch>("axisymmetricSimulation", false)),
    rWMethod_(word::null),
    revolutionAxis_(0),
    polarAxis_(1),
    angularCoordinate_(2),
    radialExtent_(0.0),
    maxRWF_(1.0),
    RWF_
    (
        IOobject
        (
            "RWF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("RWF", dimless, 1.0)
    ),
    measureMediumProperties_(particleProperties_.lookupOrDefault<Switch>("measureMediumProperties", false)),
    measureInstantaneousMSD_(particleProperties_.lookupOrDefault<Switch>("measureInstantaneousMSD", false)),
    trackingProbability_(0),
    tracerPatchNamesMap_(),
    dimensionality_(3),
    mediumTransitTime_(0.0),
    mediumTotalDistanceTravelled_(0.0),
    mediumTortuosity_(0.0),
    timeLooping_(0.0),
    nTerminalOutputs_(controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)),
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
    boundaryMeas_(mesh, *this, true),
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
    cellMeas_(mesh, *this, true)
    
{
    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();
    
    if(axisymmetric_)
    {
        rWMethod_ = particleProperties_.subDict("axisymmetricProperties")
            .lookupOrDefault<word>("radialWeightingMethod", "cell");
        
        if
        (
            rWMethod_ != "cell" and rWMethod_ != "particle"
                and rWMethod_ != "mixed"
        )
        {
            FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "The radial weighting method is badly defined. Choices "
                          "in constant/dsmcProperties are cell, particle, or "
                          "mixed. Please edit the entry: radialWeightingMethod"
                    << exit(FatalError);
        }
        
        const word& revolutionAxis = 
            particleProperties_.subDict("axisymmetricProperties")
                .lookupOrDefault<word>("revolutionAxis", word::null);
            
        const word& polarAxis = 
            particleProperties_.subDict("axisymmetricProperties")
                .lookupOrDefault<word>("polarAxis", word::null);
            
        if(revolutionAxis == "z")
        {
            revolutionAxis_ = 2;
            
            if(polarAxis == word::null) 
            {
                polarAxis_ = (revolutionAxis_ + 1)%3;
                angularCoordinate_ = (revolutionAxis_ + 2)%3;
            }
            else if(polarAxis == "y")
            {
                polarAxis_ = 1;
                angularCoordinate_ = 0;
            }
            else if(polarAxis == "x")
            {
                polarAxis_ = 0;
                angularCoordinate_ = 1;
            }
            else
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        else if(revolutionAxis == "y")
        {
            revolutionAxis_ = 1;
            
            if(polarAxis == word::null) 
            {
                polarAxis_ = (revolutionAxis_ + 1)%3;
                angularCoordinate_ = (revolutionAxis_ + 2)%3;
            }
            else if(polarAxis == "x")
            {
                polarAxis_ = 0;
                angularCoordinate_ = 2;
            }
            else if(polarAxis == "z")
            {
                polarAxis_ = 2;
                angularCoordinate_ = 0;
            }
            else
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        else if(revolutionAxis == "x")
        {
            if(polarAxis == "z")
            {
                polarAxis_ = 2;
                angularCoordinate_ = 1;
            }
            else if(polarAxis != "y")
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        
        radialExtent_ = gMax
            (
                mesh_.faceCentres().component(polarAxis_)
            );
            
        if(!(radialExtent_ > 0))
        {
            radialExtent_ = -gMin
            (
                mesh_.faceCentres().component(polarAxis_)
            );
        }
            
        maxRWF_ = readScalar
            (
                particleProperties_.subDict("axisymmetricProperties")
                    .lookup("maxRadialWeightingFactor")
            );
            
        Info << nl << "Axisymmetric simulation:" << nl
             << tab << "- revolutionAxis label" << tab << revolutionAxis_ << nl
             << tab << "- polarAxis label" << tab << polarAxis_ << nl
             << tab << "- angularCoordinate label" << tab << angularCoordinate_
             << nl << tab << "- radial weighting method" << tab << rWMethod_ 
             << "-based" << nl
             << tab << "- radialExtent" << tab << radialExtent_ << nl
             << tab << "- maxRadialWeightingFactor" << tab << maxRWF_ << nl
             << endl;
             
         updateRWF();
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
    
    dsmcParcel::TrackedParcel::nDELETED = 0;
    nLooping_ = 0;
    
    if(particleProperties_.isDict("tracerProperties"))
    {
        wordList trackFromPatchNames;
        List<wordList> trackToPatchNames;
        
        if
        (
            particleProperties_.subDict("tracerProperties")
              .found("inflowPatchNames")
        )
        {
            particleProperties_.subDict("tracerProperties")
                .lookup("inflowPatchNames") >> trackFromPatchNames;
        }
            
        if
        (
            particleProperties_.subDict("tracerProperties")
              .found("outflowPatchNames")
        )
        {
            particleProperties_.subDict("tracerProperties")
               .lookup("outflowPatchNames") >> trackToPatchNames;
        }
            
        
        if(trackFromPatchNames.size() != trackToPatchNames.size())
        {
            FatalErrorIn
            (
                "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                    "const dynamicFvMesh& mesh, bool readFields)"
            )
                << "Lists inflowPatchNames and outflowPatchNames defined"
                << " in constant/dsmcProperties must be of the same size." << nl
                << "inflowPatchNames: size " << trackFromPatchNames.size() << nl
                << "outflowPatchNames: size " << trackToPatchNames.size()
                << exit(FatalError);
        }    
            
        label noCombinations = 0;
        
        forAll(trackFromPatchNames, inflowPatch)
        {
            if(mesh_.boundaryMesh().findPatchID(trackFromPatchNames[inflowPatch]) == -1)
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "The tracers inflow patch " << trackFromPatchNames[inflowPatch]
                    << " defined in constant/dsmcProperties does not exist"
                    << exit(FatalError);
            }
            
            forAll(trackToPatchNames[inflowPatch], outflowPatch)
            {
                if(mesh_.boundaryMesh().findPatchID(trackToPatchNames[inflowPatch][outflowPatch]) == -1)
                {
                    FatalErrorIn
                    (
                        "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                            "const dynamicFvMesh& mesh, bool readFields)"
                    )
                        << "The tracers outflow patch "
                        << trackToPatchNames[inflowPatch][outflowPatch]
                        << " corresponding to the inflow patch "
                        << trackFromPatchNames[inflowPatch]
                        << " and defined in constant/dsmcProperties does not exist"
                        << exit(FatalError);
                }
                
                noCombinations++;
            }
        } 
        
        tracerPatchNamesMap_.resize(noCombinations);
        
        forAll(trackFromPatchNames, inflowPatch)
        {
            forAll(trackToPatchNames[inflowPatch], outflowPatch)
            {
                const word key = std::to_string(mesh_.boundaryMesh().findPatchID(trackFromPatchNames[inflowPatch]))
                    + std::to_string(mesh_.boundaryMesh().findPatchID(trackToPatchNames[inflowPatch][outflowPatch]));
                
                tracerPatchNamesMap_.insert(key);
            }
        }
        
        trackingProbability_ = particleProperties_
            .subDict("tracerProperties")
            .lookupOrDefault<scalar>("trackingProbability", 0.1);
            
        dimensionality_ = particleProperties_.subDict("tracerProperties")
            .lookupOrDefault<label>("dimensionality", 3);
    }
    else
    {
        if (measureMediumProperties_ or measureInstantaneousMSD_)
        {
            WarningIn
            (
                "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                    "const dynamicFvMesh& mesh, bool readFields)"
            )
                << "Subdictionary tracerProperties is missing in "
                      "constant/dsmcProperties"
                << endl;
        }
    }
    
    if(measureInstantaneousMSD_)
    {
        forAllIter(dsmcCloud, *this, iter)
        {
            if(not iter().isTracked())
            {
                if(trackingProbability_ > rndGen_.scalar01())
                {
                    iter().setTracked
                    (
                        true,
                        -1,
                        mesh_.time().value(), 
                        iter().position()
                    );
                }
            }
        }
    }
}



// running dsmcIntialise+
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
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    axisymmetric_(particleProperties_.lookupOrDefault<Switch>("axisymmetricSimulation",false)),
    rWMethod_(word::null),
    revolutionAxis_(0),
    polarAxis_(1),
    angularCoordinate_(2),
    radialExtent_(0.0),
    maxRWF_(1.0),
    RWF_
    (
        IOobject
        (
            "RWF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("RWF", dimless, 1.0)
    ),
    measureMediumProperties_(false),
    measureInstantaneousMSD_(false),
    trackingProbability_(0),
    tracerPatchNamesMap_(),
    dimensionality_(3),
    mediumTransitTime_(0.0),
    mediumTotalDistanceTravelled_(0.0),
    mediumTortuosity_(0.0),
    timeLooping_(0.0),
    nTerminalOutputs_(controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)),
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
    rndGen_(label(clock::getTime()) + 1526*Pstream::myProcNo()), // different seed every time simulation is started - needed for ensemble averaging!
    controllers_(t, mesh),
    dynamicLoadBalancing_(t, mesh, *this),
    boundaryMeas_(mesh, *this),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    cellMeas_(mesh, *this)
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
    
    if(axisymmetric_)
    {
        rWMethod_ = particleProperties_.subDict("axisymmetricProperties")
            .lookupOrDefault<word>("radialWeightingMethod", "cell");
        
        if
        (
            rWMethod_ != "cell" and rWMethod_ != "particle"
                and rWMethod_ != "mixed"
        )
        {
            FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "The radial weighting method is badly defined. Choices "
                          "in constant/dsmcProperties are cell, particle, or "
                          "mixed. Please edit the entry: radialWeightingMethod"
                    << exit(FatalError);
        }
        
        const word& revolutionAxis = 
            particleProperties_.subDict("axisymmetricProperties")
                .lookupOrDefault<word>("revolutionAxis", word::null);
            
        const word& polarAxis = 
            particleProperties_.subDict("axisymmetricProperties")
                .lookupOrDefault<word>("polarAxis", word::null);
            
        if(revolutionAxis == "z")
        {
            revolutionAxis_ = 2;
            
            if(polarAxis == word::null) 
            {
                polarAxis_ = (revolutionAxis_ + 1)%3;
                angularCoordinate_ = (revolutionAxis_ + 2)%3;
            }
            else if(polarAxis == "y")
            {
                polarAxis_ = 1;
                angularCoordinate_ = 0;
            }
            else if(polarAxis == "x")
            {
                polarAxis_ = 0;
                angularCoordinate_ = 1;
            }
            else
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        else if(revolutionAxis == "y")
        {
            revolutionAxis_ = 1;
            
            if(polarAxis == word::null) 
            {
                polarAxis_ = (revolutionAxis_ + 1)%3;
                angularCoordinate_ = (revolutionAxis_ + 2)%3;
            }
            else if(polarAxis == "x")
            {
                polarAxis_ = 0;
                angularCoordinate_ = 2;
            }
            else if(polarAxis == "z")
            {
                polarAxis_ = 2;
                angularCoordinate_ = 0;
            }
            else
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        else if(revolutionAxis == "x")
        {
            if(polarAxis == "z")
            {
                polarAxis_ = 2;
                angularCoordinate_ = 1;
            }
            else if(polarAxis != "y")
            {
                FatalErrorIn
                (
                    "dsmcCloud::dsmcCloud(Time& t, const word& cloudName, "
                        "const dynamicFvMesh& mesh, bool readFields)"
                )
                    << "Revolution and polar axes are badly defined in "
                          "constant/dsmcProperties axisymmetricProperties{}"
                    << exit(FatalError);
            }
        }
        
        radialExtent_ = gMax
            (
                mesh_.faceCentres().component(polarAxis_)
            );
            
        if(!(radialExtent_ > 0))
        {
            radialExtent_ = -gMin
            (
                mesh_.faceCentres().component(polarAxis_)
            );
        }
            
        maxRWF_ = readScalar
            (
                particleProperties_.subDict("axisymmetricProperties")
                    .lookup("maxRadialWeightingFactor")
            );
            
        Info << nl << "Axisymmetric initialisation:" << nl
             << tab << "- revolutionAxis label" << tab << revolutionAxis_ << nl
             << tab << "- polarAxis label" << tab << polarAxis_ << nl
             << tab << "- angularCoordinate label" << tab << angularCoordinate_
             << nl << tab << "- radial weighting method" << tab << rWMethod_ 
             << "-based" << nl
             << tab << "- radialExtent" << tab << radialExtent_ << nl
             << tab << "- maxRadialWeightingFactor" << tab << maxRWF_ << nl
             << endl;
             
        // "particle" cannot be used in dsmcInitialise+, "cell" is thus employed
        rWMethod_ = "cell";     
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
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();
    
    //- Remove electrons
    //removeElectrons(); // NOTE VINCENT: there is a clever way than rebuilding entire cell occ.
    //buildCellOccupancy();
    
    //- Move the particles ballistically with their current velocities
    scalar timer = mesh_.time().elapsedCpuTime(); 
    Cloud<dsmcParcel>::move(td, mesh_.time().deltaTValue());
    Info<< "move" << tab << mesh_.time().elapsedCpuTime() - timer 
        << " s" << endl; 
    
    //- Add electrons back after the move function
    //addElectrons();
    
    //- Update cell occupancy
    timer = mesh_.time().elapsedCpuTime();
    buildCellOccupancy();
    Info<< "buildCellOccupancy" << tab << mesh_.time().elapsedCpuTime() - timer
        << " s " << endl;
     
    //- Radial weighting for axially symmetric flows
    if (axisymmetric_)
    {
        if (rWMethod_ == "particle")
        {
            updateRWF();
        }
        
        axisymmetricWeighting();
        buildCellOccupancy();
    }
    
    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();
    
    //- Calculate new velocities via stochastic collisions
    timer = mesh_.time().elapsedCpuTime();
    collisions();
    Info<< "collisions" << tab << mesh_.time().elapsedCpuTime() - timer
        << " s" << endl;
    
    //- Reactions may have changed cell occupancy, update if any reaction
    if(reactions_.nReactions() != 0)
    {
        buildCellOccupancy();
    }
    
    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();
    
    reactions_.outputData();
    
    fields_.calculateFields();
    
    timer = mesh_.time().elapsedCpuTime();   
    fields_.writeFields();
    Info<< "fields W" << tab << mesh_.time().elapsedCpuTime() - timer
        << " s" << endl;

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();
    
    boundaryMeas_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean(); 
    
    //- Run-time editable entries
    measureMediumProperties_ = particleProperties_
        .lookupOrDefault<Switch>("measureMediumProperties", false);
        
    measureInstantaneousMSD_ = particleProperties_
        .lookupOrDefault<Switch>("measureInstantaneousMSD", false);    
}

Foam::label Foam::dsmcCloud::nTerminalOutputs()
{
    return nTerminalOutputs_;
}

void Foam::dsmcCloud::info()
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());
    
    const scalarList& iM = infoMeasurements();

    scalar nMol = iM[6];
    reduce(nMol, sumOp<scalar>());
    
    scalar linearKineticEnergy = iM[1];
    reduce(linearKineticEnergy, sumOp<scalar>());
    
    scalar rotationalEnergy = iM[2];
    reduce(rotationalEnergy, sumOp<scalar>());
    
    scalar vibrationalEnergy = iM[3];
    reduce(vibrationalEnergy, sumOp<scalar>());
    
    scalar electronicEnergy = iM[4];
    reduce(electronicEnergy, sumOp<scalar>());
    
    scalar stuckMolecules = iM[5];
    reduce(stuckMolecules, sumOp<scalar>());

    Info<< "    Number of DSMC particles        = "
        << nDsmcParticles
        << endl;

    if (nDsmcParticles > VSMALL)
    {
       Info << "    Number of stuck particles       = "
            << stuckMolecules/nParticle_ << nl
            << "    Number of free particles        = "
            << nMol/nParticle_ << nl
            << "    Average linear kinetic energy   = "
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
            
            if (tracerPatchNamesMap_.size() != 0)
            {
                const scalar nDel = 
                    dsmcParcel::TrackedParcel::nDELETED + SMALL;
                const scalar nLooping = nLooping_ + SMALL;
                
                scalar weightingFactor = 0.0;
                
                if(dsmcParcel::TrackedParcel::nDELETED > 0)
                {
                    weightingFactor = nLooping_/nDel;
                }
                
                const scalar mediumSpatialExtension = 1.0;
                
                mediumTortuosity_ = mediumTotalDistanceTravelled_
                    /(nDel*mediumSpatialExtension);
                
                Info << "    Tracked parcels recorded        = " 
                     << dsmcParcel::TrackedParcel::nDELETED << nl
                     << "    Mean particle displacement      = " 
                     << mediumTotalDistanceTravelled_/nDel << nl
                     << "    Medium tortuosity               = " 
                     << mediumTortuosity_ << nl
                     << "    Recorded flow transit time      = "
                     << mediumTransitTime_/nDel << nl
                     << "    Looping parcels recorded        = " 
                     << nLooping_ << nl
                     << "    Flow time penalty               = "
                     << timeLooping_/nLooping << nl 
                     << "    Corrected flow transit time     = "
                     << mediumTransitTime_/nDel 
                        + weightingFactor*timeLooping_/nLooping
                     << endl;
             }
             
             if (measureInstantaneousMSD_)
             {
                  vector iMSDvector = vector::zero;
                  label cmpt = 0;
                  
                  forAllIter(dsmcCloud, *this, iter)
                  {
                      dsmcParcel& p = iter();
                      
                      if (p.isTracked())
                      {
                          if (p.tracked().inPatchId() == -1)
                          {
                              p.tracked().updateDistanceTravelled(p.position());
                              
                              p.tracked().updateMeanSquaredDisplacement();
                              
                              iMSDvector += 
                                  p.tracked().meanSquaredDisplacementVector();
                              
                              cmpt++;
                          }
                      }
                  }
                  
                  scalar iMSD = iMSDvector[0] + iMSDvector[1] + iMSDvector[2];
                    
                  Info<< "    Instantaneous msd               = " 
                      << iMSD/(cmpt + SMALL) << " (" << cmpt << ")" << nl
                      << "    IMSDx                           = " 
                      << iMSDvector[0]/(cmpt + SMALL) << nl
                      << "    IMSDy                           = " 
                      << iMSDvector[1]/(cmpt + SMALL) << nl
                      << "    IMSDz                           = " 
                      << iMSDvector[2]/(cmpt + SMALL)
                      << endl;
             }
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

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();
    resetBoundaries();
    resetMeasurementTools();
}


Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar& temperature,
    const scalar& mass
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)
       *vector
        (
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal(),
            rndGen_.GaussNormal()
        );
}


Foam::vector Foam::dsmcCloud::chapmanEnskogVelocity
(
    const scalar& temperature,
    const scalar& mass,
    const vector& q,
    const tensor& tau
)
{
    const scalar B = max(mag(q), mag(tau));
    const scalar A = 1.0 + 30.0*B;

    bool repeatTry = true;
    
    vector CTry = vector::zero;
    
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
            
        if (gammaTry >= A*rndGen_.scalar01())
        {
            repeatTry = false;
        }
    }
    
    return CTry*sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


void Foam::dsmcCloud::generalisedChapmanEnskog
(
    const label& typeID,
    const scalar& translationalTemperature,
    const scalar& rotationalTemperature,
    const scalar& vibrationalTemperature,
    const scalar& mass,
    const vector& D,
    const vector& qTra,
    const vector& qRot,
    const vector& qVib,
    const tensor& tau,
    scalar& ERot,
    labelList& vibLevel,
    vector& U
)
{
    const scalar k = physicoChemical::k.value();
    const scalar vibDOF = constProps(typeID).vibrationalDegreesOfFreedom();

    scalar B = max(mag(D), mag(tau));
    B = max(B, mag(qTra));
    B = max(B, mag(qRot));
    B = max(B, mag(qVib));
    const scalar A = 1.0 + 30.0*B;

    const scalarList& thetaV = constProps(typeID).thetaV();
    const scalar epsRotAv = rotationalTemperature/translationalTemperature;
    
    scalar epsVibAv = 0.0;
    
    if(vibDOF > 0 && vibrationalTemperature > SMALL)
    {
        forAll(vibLevel, i)
        {
            epsVibAv += thetaV[i]/vibrationalTemperature
                /(exp(thetaV[i]/vibrationalTemperature) - 1.0);
        }
        
    }

    vector CTry = vector::zero;
    
    scalar epsRot = 0.0;
    scalar epsVib = 0.0;
    
    bool repeatTry = true;
    
    while (repeatTry)
    {
        ERot = equipartitionRotationalEnergy
            (
                rotationalTemperature,
                constProps(typeID).rotationalDegreesOfFreedom()
            );
        
        vibLevel = equipartitionVibrationalEnergyLevel
            (
                vibrationalTemperature,
                vibDOF,
                typeID
            );


        epsRot = ERot/(k*translationalTemperature);
        
        if(vibDOF > 0 && vibrationalTemperature > SMALL)
        {
            forAll(vibLevel, i)
            {
                epsVib += vibLevel[i]
                    *physicoChemical::k.value()
                    *thetaV[i];
            }
            
            epsVib /= (physicoChemical::k.value()*vibrationalTemperature);
        }

        CTry = vector
            (
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal(),
                rndGen_.GaussNormal()
            ) / sqrt(2.0);
            
        scalar gammaTry = 1.0 + 2.0*(D & CTry) 
            + (qTra & CTry) * (0.4*(CTry & CTry) - 1.0) 
            + (qRot & CTry) * (epsRot - epsRotAv)
            + (qVib & CTry) * (epsVib - epsVibAv)
            - (CTry & (tau & CTry)); 
            
        if (gammaTry >= A*rndGen_.scalar01())
        {
            repeatTry = false;
        }
    }

    U = CTry*sqrt(2.0*k*translationalTemperature/mass);
}


Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    const scalar& temperature,
    const scalar& rotationalDof
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
    const scalar& temperature,
    const scalar& vibrationalDof,
    const label& typeId
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
    const scalar& temperature,
    const List<label>& degeneracyList,
    const List<scalar>& electronicEnergyList,
    const label& typeId
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
            expSum += degeneracyList[i]*exp((-electronicEnergyList[i]/EMax));  
            i += 1;
        } while (i < jMax);     

        // select maximum integer energy level based on boltz value.
        // Note that this depends on the temperature.   
        
        scalar boltzMax = 0.0;
        
        for (label ii = 0; ii < jMax; ii++)
        {
            //Eq. 3.1.1 of Liechty thesis.   
            boltz = degeneracyList[ii]*exp((-electronicEnergyList[ii]/EMax))/expSum;
            
            if (boltzMax < boltz)
            {
                boltzMax = boltz;
                jSelect = ii;
            }               
        }

        EJ = electronicEnergyList[jSelect]; //Max. poss energy in list : list goes from [0] to [jMax-1]
        gJ = degeneracyList[jSelect]; //Max. poss degeneracy in list : list goes from [0] to [jMax-1]
        expMax = gJ*exp((-EJ/EMax)); // Max. in denominator of Liechty pdf for initialisation/
                                     //wall bcs/freestream EEle etc..
        
        do // acceptance - rejection based on Eq. 3.1.2 of Liechty thesis.      
        {               
            jDash = rndGen_.integer(0,jMax-1);
            func = degeneracyList[jDash]*exp((-electronicEnergyList[jDash]/EMax))/expMax;
        } while( !(func > rndGen_.scalar01()));
    }

    return jDash;
}


Foam::scalar Foam::dsmcCloud::postCollisionRotationalEnergy
(
    const scalar& rotationalDof,
    const scalar& ChiB  
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
    scalar Ec,
    scalar fixedZv
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

        scalar inverseVibrationalCollisionNumber = 1.0; // NEW VINCENT
        
        if (fixedZv == 0)
        {
            //- vibrational collision number (equation 2, Bird 2010)
            scalar ZvP1 = pow((thetaD/TColl),omega); 
        
            scalar ZvP2 = pow(Zref*(pow((thetaD/refTempZv),(-1.0*omega))),(pow1/pow2));
            
            scalar Zv = ZvP1*ZvP2;
            
            //- In order to obtain the relaxation rate corresponding to Zv with the collision 
            //  energy-based procedure, the inelastic fraction should be set to about 1/(5Zv)
            //  Bird 2008 RGD "A Comparison of Collision Energy-Based and Temperature-Based..."
            
            inverseVibrationalCollisionNumber = 1.0/(5.0*Zv);
        }
        else
        {
            inverseVibrationalCollisionNumber = 1.0/fixedZv;
        }
        
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
    forAll(cellOccupancy_, c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            
            const scalar oldRadialWeight = p->RWF();
                        
            const scalar newRadialWeight = RWF(c); // getRWF_cell(c);

            p->RWF() = newRadialWeight;
            
            if(oldRadialWeight > newRadialWeight) 
            {
                //- particle might be cloned
                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;
                
                while(prob > 1.0)
                {
                    //- add a particle and reduce prob by 1.0
                    vector U = p->U();
                    
                    //U.z() *= -1.0;
                    U.component(angularCoordinate_) *= -1.0;
                    
                    addNewParcel
                    (
                        p->position(),
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        p->cell(),
                        p->tetFace(),
                        p->tetPt(),
                        p->typeId(),
                        p->newParcel(),
                        p->classification(),
                        p->vibLevel()
                    );
                    
                    prob -= 1.0;
                }
                
                if(prob > rndGen_.scalar01())
                {
                    vector U = p->U();
                    
                    //U.z() *= -1.0;
                    U.component(angularCoordinate_) *= -1.0;

                    addNewParcel
                    (
                        p->position(),
                        U,
                        p->RWF(),
                        p->ERot(),
                        p->ELevel(),
                        p->cell(),
                        p->tetFace(),
                        p->tetPt(),
                        p->typeId(),
                        p->newParcel(),
                        p->classification(),
                        p->vibLevel()
                    );
                }
            }
            else if(newRadialWeight > oldRadialWeight)
            {           
                //- particle might be deleted
                if((oldRadialWeight/newRadialWeight) < rndGen_.scalar01())
                {
                    deleteParticle(*p);
                } 
            } 
        }
    }
}


void Foam::dsmcCloud::updateRWF()
{
    forAll(RWF_, celli)
    {
        RWF_[celli] = getRWF_cell(celli);
    }
    
    forAll(RWF_.boundaryField(), patchi)
    {
        fvPatchScalarField& pRWF = RWF_.boundaryFieldRef()[patchi];
        
        forAll(pRWF, facei)
        {
            pRWF[facei] = getRWF_face(patchi, facei);
        }
    }
}


Foam::scalar Foam::dsmcCloud::getRWF_face
(
    const label patchI,
    const label faceI
) const
{
    scalar RWF = 1.0;
    
    if(axisymmetric_)
    {
        const point& fC = mesh_.boundaryMesh()[patchI].faceCentres()[faceI];
        const scalar radius = mag(fC.component(polarAxis_));
        
        RWF += (maxRWF() - 1.0)*radius/radialExtent();
    }

    return RWF;    
}


Foam::scalar Foam::dsmcCloud::getRWF_cell
(
    const label cellI, 
    const bool mixedRWM
) const
{
    scalar RWF = 1.0;
    
    if(axisymmetric_)
    {
        if(rWMethod_ == "particle" or (mixedRWM and rWMethod_ == "mixed"))
        {
            const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy_[cellI]);
            
            RWF = 0.0;
            label nMols = 0;
            
            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];
                
                const scalar radius = 
                    sqrt
                    (
                        sqr(p.position().component(polarAxis_))
                      + sqr(p.position().component(angularCoordinate_))
                    );

                RWF += 1.0 + (maxRWF() - 1.0)*radius/radialExtent();
                
                nMols++;
            }
            
            RWF /= max(nMols, 1);
        }
        else
        {
            const point& cC = mesh_.cellCentres()[cellI];
            const scalar radius = mag(cC.component(polarAxis_));
        
            RWF += (maxRWF() - 1.0)*radius/radialExtent();
        }
    }
    
    return RWF;    
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


// NEW DANIEL *****************************************************************
/*void Foam::dsmcCloud::resetHybrid
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
                
                            label newParcel = -1;
            
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
                
                            label newParcel = -1;
            
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
                
                            label newParcel = -1;
            
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
                
                            label newParcel = -1;
            
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
}*/

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
    molsToDelete(mesh_, *this, typeOfReset);
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

                            scalar ERot = 0.0;
                            
                            labelList vibLevel
                            (
                                constProps(i).thetaV().size(),
                                0.0
                            );
                            
                            label ELevel = 0; // TODO by generalisedChapmanEnskog

                            generalisedChapmanEnskog
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                constProps(i).mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                vibLevel,
                                U
                            );

                            U += velocity;
                
                            label newParcel = -1;
            
                            label classification = 0;
                            
                            const scalar& RWF = RWF_[cellI]; //getRWF_cell(cellI);

                            addNewParcel
                            (
                                p,
                                U,
                                RWF,
                                ERot,
                                ELevel,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification,
                                vibLevel
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

/*void Foam::dsmcCloud::resetHybridTraRotVib2
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
                
                            label newParcel = -1;
            
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
}*/


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
// END NEW DANIEL *************************************************************


void Foam::dsmcCloud::updateMediumPropertiesMeasurement
(
    dsmcParcel& p,
    const label& delPatchId
)
{
    if(measureMediumProperties_)
    {
        if(p.isTracked())
        {
            if (p.tracked().inPatchId() != -1)
            {
                p.tracked().updateTotalDistanceTravelled(p.position());
                
                const word patchesKey = std::to_string(p.tracked().inPatchId()) 
                    + std::to_string(delPatchId);
                
                if(tracerPatchNamesMap_.found(patchesKey))
                {
                    //- The deletion patch is in the list of outflow patches of the
                    //  tracer inflow patch
                    dsmcParcel::TrackedParcel::nDELETED++;
                    
                    mediumTotalDistanceTravelled_ += 
                        p.tracked().distanceTravelled();
                    
                    mediumTransitTime_ += mesh_.time().value()
                        - p.tracked().initialTime();
                }
                else if(patchesKey[0] == patchesKey[1])
                {
                    //- The deletion patch is NOT in the list of outflow patches of
                    //  the tracer inflow patch. Add a time penalty
                    nLooping_++;
                    
                    timeLooping_ += mesh_.time().value() 
                        - p.tracked().initialTime();
                }
            }
        }
    }
}


void Foam::dsmcCloud::updateMediumPropertiesMeasurement_cyclic
(
    dsmcParcel& p,
    const label& neiPatchId,
    const label& orgPatchId,
    const vector& orgPosition
)
{
    /*if(measureMediumProperties_)
    {
        if(p.isTracked())
        {
            // orgPatchId is the id of the patch on which the particle strikes, i.e., the 'deletion' patch
            // neiPatchId is the id of the neighbouring patch, i.e., the patch of reintroduction
            
            // TODO update mediumTotalDistanceTravelled_
            
            if(orgPatchId != p.tracked().inPatchId())
            {
                //- The cyclic patch is different from the patch of insertion
                //Info << "inP: " << p.tracked().inPatchId() << tab << "orgP: " << orgPatchId << tab << "neiP: " << neiPatchId << endl;
                dsmcParcel::TrackedParcel::nDELETED++;
                
                mediumTotalDistanceTravelled_ += p.tracked().distanceTravelled();
                    
                mediumTransitTime_ += mesh_.time().value() - p.tracked().initialTime();
                
                if(trackingProbability_ > rndGen_.scalar01())
                {
                    //- Since the DSMC parcel is not physically deleted,
                    //  the TrackedParcel is reset.
                    p.tracked().setInitialParcelInfo
                    (
                        neiPatchId,
                        mesh_.time().value(),
                        p.position()
                    );
                }
                else
                {
                    //- Otherwise, the tracked data is deleted
                    p.deleteTracked();
                }
            }
            else
            {
                //- The patch of deletion is identical to the patch of insertion
                //  1) Check if this cyclic patch is present in the list of inflow patches
                label patchNameCounter = 0;
                do
                {
                    if(orgPatchId == mesh_.boundaryMesh().findPatchID(trackFromPatchNames_[patchNameCounter]))
                    {
                        // 1.1) If so, add a time penalty
                        nLooping_++;
                        timeLooping_ += mesh_.time().value() - p.tracked().initialTime();
                        break;
                    }
                    
                    patchNameCounter++;
                    
                } while(patchNameCounter < trackFromPatchNames_.size());
                
                // 2) Delete the tracked data
                p.deleteTracked();
            }
        }
        else
        {
            //- Every particle crossing a cyclic patch has a chance to be 
            //  tracked since there is no inlet for that purpose.
            if(trackingProbability_ > rndGen_.scalar01())
            {   
                p.setTracked
                (
                    true, 
                    neiPatchId,
                    mesh_.time().value(),
                    p.position()
                );
            }
        }
    }*/
}


bool Foam::dsmcCloud::read()
{
    if(regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
