/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "polyMoleculeCloud.H"
#include "polyAllConfigurations.H"
#include "fvMesh.H"
#include "polyMolsToDelete.H"
#include "polyMappingModels.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<polyMolecule>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
/*
void Foam::polyMoleculeCloud::buildConstProps()
{
    Info<< nl << "Reading moleculeProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size()); 

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    dictionary moleculeProperties
    (
        moleculePropertiesDict.subDict("moleculeProperties")
    );

    forAll(idList, i)
    {
        const word& id(idList[i]);

        const dictionary& molDict(moleculeProperties.subDict(id));

        const word cloudType = molDict.lookup("cloudType");

        if(cloudType == "polyMoleculeCloud")
        {
            List<word> siteIdNames = molDict.lookup("siteIds");
    
            List<label> siteIds(siteIdNames.size());
    
            forAll(siteIdNames, sI)
            {
                const word& siteId = siteIdNames[sI];
    
                siteIds[sI] = findIndex(siteIdList, siteId);
    
                if (siteIds[sI] == -1)
                {
                    FatalErrorIn("polyMoleculeCloud.C") << nl
                        << siteId << " site not found."
                        << nl << abort(FatalError);
                }
            }
    
            polyMolecule::constantProperties& constProp = constPropList_[i];
    
            constProp = polyMolecule::constantProperties(molDict, redUnits_, siteIds);
        }
    }
}
*/

void Foam::polyMoleculeCloud::setSiteSizesAndPositions()
{
    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
//         const polyMolecule::constantProperties& cP = constProps(mol().id());
        
        mol().setSiteSizes(cP_.nSites(mol().id()));

        mol().setSitePositions(cP_);
    }
}

void Foam::polyMoleculeCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    iterator mol(this->begin());

    for
    (
        mol = this->begin();
        mol != this->end();
        ++mol
    )
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }
    
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }
}

// NEW //
void Foam::polyMoleculeCloud::checkForOverlaps()
{
    if(p_.checkPotentialOverlaps())
    {
        const scalar& potLim = p_.potentialEnergyLimit();
        
        Info<< nl << "Removing high energy overlaps, limit = "
            << potLim
            << ", from removalOrder list = " << p_.removalOrder()
            << endl;

        label initialSize = this->size();

        if (Pstream::parRun())
        {
            reduce(initialSize, sumOp<label>());
        }

        buildCellOccupancy();

        prepareInteractions();

        DynamicList<polyMolecule*> molsToDelete;
        DynamicList<label> molsToDeleteTNs;

        polyMolecule* molI = NULL;
        polyMolecule* molJ = NULL;


        {
            // Real-Real interactions
            const labelListList& dil = iL_.dil();

            forAll(dil, d)
            {
                forAll(cellOccupancy_[d],cellIMols)
                {
                    molI = cellOccupancy_[d][cellIMols];
                    label idI = molI->id();
                    bool molIDeleted = false;
                    label tNI=molI->trackingNumber();

                    forAll(dil[d], interactingCells)
                    {
                        List<polyMolecule*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                        forAll(cellJ, cellJMols)
                        {
                            molJ = cellJ[cellJMols];
                            label tNJ = molJ->trackingNumber();

                            label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                            if(!molIDeleted && (molJDeleted == -1))
                            {
                                if(evaluatePotentialLimit(molI, molJ, potLim))
                                {
                                    label idJ = molJ->id();

                                    label removeIdI = findIndex(p_.removalOrder(), idI);
                                    label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                    if(removeIdI < removeIdJ)
                                    {
                                        molsToDelete.append(molI);
                                        molsToDeleteTNs.append(tNI);
                                        molIDeleted = true;
                                    }
                                    else
                                    {
                                        molsToDelete.append(molJ);
                                        molsToDeleteTNs.append(tNJ);
                                    }
                                }
                            }
                        }
                    }

                    forAll(cellOccupancy_[d], cellIOtherMols)
                    {
                        molJ = cellOccupancy_[d][cellIOtherMols];
                        label tNJ = molJ->trackingNumber();

                        label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                        if(!molIDeleted && (molJDeleted == -1))
                        {
                            if (molJ > molI)
                            {
                                if(evaluatePotentialLimit(molI, molJ, potLim))
                                {
                                    label idJ = molJ->id();

                                    label removeIdI = findIndex(p_.removalOrder(), idI);
                                    label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                    if(removeIdI < removeIdJ)
                                    {
                                        molsToDelete.append(molI);
                                        molsToDeleteTNs.append(tNI);
                                        molIDeleted = true;
                                    }
                                    else
                                    {
                                        molsToDelete.append(molJ);
                                        molsToDeleteTNs.append(tNJ);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        {
            // Real-Referred interactions
            forAll(iL_.refCellsParticles(), r)
            {
                const List<label>& realCells = iL_.refCells()[r].neighbouringCells();

                forAll(iL_.refCellsParticles()[r], i)
                {
                    molJ = iL_.refCellsParticles()[r][i];
                    label tNJ = molJ->trackingNumber();

                    label molJDeleted = findIndex(molsToDeleteTNs, tNJ);

                    if(molJDeleted == -1)
                    {
                        forAll(realCells, rC)
                        {
                            List<polyMolecule*> molsInCell = cellOccupancy_[realCells[rC]];

                            forAll(molsInCell, j)
                            {
                                molI = molsInCell[j];
                                label tNI = molI->trackingNumber();

                                label molIDeleted = findIndex(molsToDeleteTNs, tNI);

                                if(molIDeleted == -1)
                                {
                                    if(evaluatePotentialLimit(molI, molJ, potLim))
                                    {
                                        label idJ = molJ->id();
                                        label idI = molI->id();

                                        label removeIdI = findIndex(p_.removalOrder(), idI);
                                        label removeIdJ = findIndex(p_.removalOrder(), idJ);

                                        if(removeIdI < removeIdJ)
                                        {
                                            molsToDelete.append(molI);
                                            molsToDeleteTNs.append(tNI);
                                        }
                                        else if(removeIdI == removeIdJ)
                                        {
                                            if (molI->trackingNumber() > molJ->trackingNumber())
                                            {
                                                molsToDelete.append(molI);
                                                molsToDeleteTNs.append(tNI);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        label nMolsDeleted = 0;

        forAll (molsToDelete, mTD)
        {
            nMolsDeleted++;

            Pout << nl << " WARNING: Deleting molecule "
                <<  " proc no = " << Pstream::myProcNo()
                << ", position = " << molsToDelete[mTD]->position()
                << ", molecule = " << cP_.molIds()[molsToDelete[mTD]->id()]
                << endl;

            deleteParticle(*(molsToDelete[mTD]));

        }

        if (Pstream::parRun())
        {
            reduce(nMolsDeleted, sumOp<label>());
        }

        if(nMolsDeleted > 0)
        {
            Info << nl << " WARNING: Total number of molecules deleted = " << nMolsDeleted << endl;
        }
        else
        {
            Info << " NO OVERLAPPING MOLECULES" << endl;
        }

        molsToDelete.clear();
    }
    
    buildCellOccupancy();

    prepareInteractions();
}
/*
void Foam::polyMoleculeCloud::removeHighEnergyOverlaps()
{
    Info<< nl << "Removing high energy overlaps, limit = "
        << pot_.potentialEnergyLimit()
        << nl << "Removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        if(pot_.removalOrder()[rO] != -1)
        {
            Info<< ' ' << pot_.idList()[pot_.removalOrder()[rO]];
        }
    }

    Info<< nl ;

    label initialSize = this->size();

    if (Pstream::parRun())
    {
        reduce(initialSize, sumOp<label>());
    }

    buildCellOccupancy();
    
    label nMolsDeleted = 0;

    prepareInteractions();

    setIPL();
    
    iL_.setRIPL();
    
    label nMolsInt = 0;
    label nMolsExt = 0;
    label nMolsRef = 0;

    {
        DynamicList<polyMolecule*> molsToDelete;

        polyMolecule* molI = NULL;
        polyMolecule* molJ = NULL;

        forAll(ipl_, c)
        {
            nMolsExt = ipl_[c].size();
            nMolsInt = cellOccupancy_[c].size();
            nMolsRef = iL_.ripl()[c].size();
                 
            for (int i = 0; i < nMolsInt; i++)
            {
                molI = cellOccupancy_[c][i];

				label idI = molI->id();
				bool molIDeleted = false;

				for (int j = 0; j < nMolsInt; j++)
				{
					if(j > i)
					{
						molJ = cellOccupancy_[c][j];

						label molJDeleted = findIndex(molsToDelete, molJ);

						if(!molIDeleted && (molJDeleted == -1))
						{
							if
							(
								evaluatePotentialLimit
								(
									molI,
									molJ,
									pot_.potentialEnergyLimit()
								)
							)
							{
								label idJ = molJ->id();

								label removeIdI = findIndex(pot_.removalOrder(), idI);
								label removeIdJ = findIndex(pot_.removalOrder(), idJ);

								if(removeIdI < removeIdJ)
								{
									molsToDelete.append(molI);
									molIDeleted = true;
								}
								else
								{
									molsToDelete.append(molJ);
								}
							}
						}
					}
				}

				for (int j = 0; j < nMolsExt; j++)
				{
					molJ = ipl_[c][j];

					label molJDeleted = findIndex(molsToDelete, molJ);

					if(!molIDeleted && (molJDeleted == -1))
					{
						if
						(
							evaluatePotentialLimit
							(
								molI,
								molJ,
								pot_.potentialEnergyLimit()
							)
						)
						{
							label idJ = molJ->id();

							label removeIdI = findIndex(pot_.removalOrder(), idI);
							label removeIdJ = findIndex(pot_.removalOrder(), idJ);

							if(removeIdI < removeIdJ)
							{
								molsToDelete.append(molI);
								molIDeleted = true;
							}
							else
							{
								molsToDelete.append(molJ);
							}
						}
					}
				}

				for (int j = 0; j < nMolsRef; j++)
				{
					molJ = iL_.ripl()[c][j];

					if(!molIDeleted)
					{
						if
						(
							evaluatePotentialLimit
							(
								molI,
								molJ,
								pot_.potentialEnergyLimit()
							)
						)
						{
							label idJ = molJ->id();

							label removeIdI = findIndex(pot_.removalOrder(), idI);
							label removeIdJ = findIndex(pot_.removalOrder(), idJ);

							if(removeIdI < removeIdJ)
							{
								molsToDelete.append(molI);
								molIDeleted = true;
							}
							else if
							(
								(removeIdI == removeIdJ)
							)
							{
								if (molI->trackingNumber() > molJ->trackingNumber())
								{
									molsToDelete.append(molI);
									molIDeleted = true;
								}
							}
						}
					}
				}
            }
        }

        forAll (molsToDelete, mTD)
        {
            nMolsDeleted++;
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    label newSize = this->size();

    if (Pstream::parRun())
    {
        reduce(newSize, sumOp<label>());
    }
    
    if (Pstream::parRun())
    {
        reduce(nMolsDeleted, sumOp<label>());
    }    

    if(nMolsDeleted > 0)
    {
        // to make sure the user sees this
        for (int j = 0; j < 50; j++)
        {
            Info << nl << "WARNING: molecules removed due to overlaps = "
                << nMolsDeleted <<  endl;
        }
    }
    
}*/


Foam::label Foam::polyMoleculeCloud::nSites() const
{
    label n = 0;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
//         n += constProps(mol().id()).nSites();
        n += cP_.nSites(mol().id());
    }

    return n;
}


void Foam::polyMoleculeCloud::checkMoleculesInMesh()
{
    Info << nl << "checking cell-molecule addressing" << endl;

    DynamicList<polyMolecule*> molsToDelete;

    label initialSize = this->size();

    iterator mol(this->begin());

    label noOfModifiedMols = 0;

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            mol().position(),
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            if(mol().cell() != cell)
            {
                mol().cell() = cell;
                mol().tetFace() = tetFace;
                mol().tetPt() = tetPt;
                noOfModifiedMols++;
            }
        }
        else
        {
            polyMolecule* molI = &mol();
            molsToDelete.append(molI);
        }
    }

    if(noOfModifiedMols > 0)
    {
        Pout<< tab << " molecules that changed cell = " 
            << noOfModifiedMols
            << endl;
    }

    forAll (molsToDelete, mTD)
    {
        Pout << nl << " WARNING: Molecule Outside Mesh - Deleting molecule "
                    <<  " proc no = " << Pstream::myProcNo()
                    << ", position = " << molsToDelete[mTD]->position()
                    << ", molecule = " << cP_.molIds()[molsToDelete[mTD]->id()]
                    << endl;

        deleteParticle(*(molsToDelete[mTD]));
    }

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab <<" molecules removed from outside mesh = " 
        << molsRemoved 
        << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


//- Use for running MD (mdFoam)
Foam::polyMoleculeCloud::polyMoleculeCloud
(
    Time& t,
    const polyMesh& mesh,
//     const potentials& p,
    const reducedUnits& rU,
    const constantMoleculeProperties& cP, 
    cachedRandomMD& rndGen
)
:
    Cloud<polyMolecule>(mesh, "polyMoleculeCloud", false),
    mesh_(mesh),
//     p_(p),
    redUnits_(rU),
    cP_(cP),
    rndGen_(rndGen),
    int_(t, mesh_, *this),    
    p_(mesh, *this, rU, cP), 
    cellOccupancy_(mesh_.nCells()),
//     constPropList_(),
    fields_(t, mesh_, *this),
    boundaries_(t, mesh, *this),
    controllers_(t, mesh, *this),
    trackingInfo_(mesh, *this),
    moleculeTracking_(),
    cyclics_(t, mesh_, -1), 
    iL_(mesh, rU, cyclics_, p_.rCutMax(), "poly"),
    ipl_(mesh.nCells()),
	clock_(t, "evolve", true)
{
    polyMolecule::readFields(*this);

    rndGen.initialise(this->size() != 0 ? this->size() : 10000); //Initialise the random number cache (initialise to 10000 if size is zero)

//     buildConstProps();
    
    setSiteSizesAndPositions();

    checkMoleculesInMesh();

    // read in tracking numbers
    updateTrackingNumbersAfterRead();
    p_.pairPots().initialiseExclusionModels();

    int_.integrator()->init();
    
    //check and remove high energy overalps
    checkForOverlaps();
//     controllers_.initialConfig();
    
    buildCellOccupancy();

    
    fields_.createFields();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
    
    clearLagrangianFields();
    calculateForce();
    updateAcceleration();
    

    
    // TESTS
    writeReferredCloud();
}



//- general constructor
Foam::polyMoleculeCloud::polyMoleculeCloud
(
    Time& t,
    const polyMesh& mesh,
//     const potentials& p,
    const reducedUnits& rU,
    const constantMoleculeProperties& cP,
    cachedRandomMD& rndGen, 
    const word& option,
    const bool& clearFields
)
    :
    Cloud<polyMolecule>(mesh, "polyMoleculeCloud", false),
    mesh_(mesh),
    redUnits_(rU),
    cP_(cP),
    rndGen_(rndGen),    
    int_(t, mesh_, *this),    
    p_(mesh, *this, rU, cP), 
    cellOccupancy_(mesh_.nCells()),
//     constPropList_(),
    fields_(t, mesh_),
    boundaries_(t, mesh),
    controllers_(t, mesh),
    trackingInfo_(mesh, *this),
    moleculeTracking_(),
    cyclics_(t, mesh_, -1),
    iL_(mesh, rU, cyclics_, p_.rCutMax(), "poly"),
    ipl_(mesh.nCells()),
	clock_(t, "evolve", true)
{
    polyMolecule::readFields(*this);

    label initialMolecules = this->size();

    rndGen.initialise(initialMolecules != 0 ? initialMolecules : 10000); //Initialise the random number cache (initialise to 10000 if size is zero)

    if (Pstream::parRun())
    {
        reduce(initialMolecules, sumOp<label>());
    }
   
    if(clearFields)
    {
        Info << "clearing existing field of molecules " << endl;

        clear();

        initialMolecules = 0;
    }
    else
    {
        updateTrackingNumbersAfterRead();
        setSiteSizesAndPositions();        
    }

    if((option == "mdInitialise") && clearFields)
    {
        polyAllConfigurations conf(mesh, *this);
        conf.setInitialConfig();
        buildCellOccupancy();
    }
    else if((option == "mdInitialise") && !clearFields)
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
        polyAllConfigurations conf(mesh, *this);
        conf.setInitialConfig();
    }
    else if(option == "delete")
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
        prepareInteractions();
        polyMolsToDelete molsDel(mesh_, *this);
    }
    else if(option == "mapping")
    {
        polyMappingModels molsToMap(mesh_, *this);
        buildCellOccupancy();
    }
    else if(option == "quickMapping")
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
    }
    else if(option == "NULL")
    {
        buildCellOccupancy();
    }
    else 
    {
        Info << "ERROR" << endl;
    }

    label finalMolecules = this->size();
    
    if (Pstream::parRun())
    {
        reduce(finalMolecules, sumOp<label>());
    }

    Info << nl << "Initial molecules = " << initialMolecules 
         << ", modified molecules = " << finalMolecules - initialMolecules
         << ", total molecules: " << finalMolecules 
         << endl;
}

// * * * * * * * * * * * * * * * * Static Constructors  * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::polyMoleculeCloud> Foam::polyMoleculeCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const constantMoleculeProperties& cP, 
    cachedRandomMD& rndGen
)
{
    return autoPtr<polyMoleculeCloud>
    (
        new polyMoleculeCloud(t, mesh, rU, cP, rndGen)
    );
}

Foam::autoPtr<Foam::polyMoleculeCloud> Foam::polyMoleculeCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const reducedUnits& rU,
    const constantMoleculeProperties& cP, 
    cachedRandomMD& rndGen,
    const word& option,
    const bool& clearFields
)
{
    return autoPtr<polyMoleculeCloud>
    (
        new polyMoleculeCloud(t, mesh, rU, cP, rndGen, option, clearFields)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void  Foam::polyMoleculeCloud::createMolecule
(
    const vector& position,
    const label cell,
    const label tetFace,
    const label tetPt,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    const vector& specialPosition,
    const label special,
    const label id,
    const scalar& fraction,
    const label trackingNumber
)
{
    addParticle
    (
        new polyMolecule
        (
            mesh_,
            position,
            cell,
            tetFace,
            tetPt,
            Q,
            v,
            a,
            pi,
            tau,
            specialPosition,
//             constProps(id),
            cP_,
            special,
            id,
            fraction,
            trackingNumber
        )
    );
}


// Evolve functions 

void Foam::polyMoleculeCloud::evolve()
{
    int_.integrator()->evolve();

//     evolveBeforeForces();
//     calculateForce();
//     evolveAfterForces();
}

void Foam::polyMoleculeCloud::evolveBeforeForces()
{
    controlBeforeVelocity();
    updateVelocity();
    controlBeforeMove();
    move();
    controlAfterMove();
    buildCellOccupancy();
    controlBeforeForces();
    clearLagrangianFields();
}

void Foam::polyMoleculeCloud::evolveAfterForces()
{
    updateAcceleration();
    controlAfterForces();
    updateVelocity();
    controlAfterVelocity();
    postTimeStep();
}

void Foam::polyMoleculeCloud::controlBeforeVelocity()
{
    controllers_.controlVelocitiesI();
}

void Foam::polyMoleculeCloud::updateVelocity()
{
    velocityUpdate(mesh_.time().deltaT().value());
}

void Foam::polyMoleculeCloud::velocityUpdate(const scalar& trackTime)
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            mol().updateHalfVelocity(cP_, trackTime);
        }
    }
}

void Foam::polyMoleculeCloud::controlBeforeMove()
{
    controllers_.controlBeforeMove();
}

// move molecules (tracking)
void Foam::polyMoleculeCloud::move()
{
    polyMolecule::trackingData td1(*this, 1);
    Cloud<polyMolecule>::move(td1, mesh_.time().deltaTValue());

    updateAfterMove(mesh_.time().deltaT().value());
}

void Foam::polyMoleculeCloud::move(const scalar& trackTime)
{
    polyMolecule::trackingData td1(*this, 1);
    Cloud<polyMolecule>::move(td1, trackTime);
}


void Foam::polyMoleculeCloud::updateAfterMove(const scalar& trackTime)
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            mol().updateAfterMove(cP_, trackTime);
        }
    }
}


void Foam::polyMoleculeCloud::controlAfterMove()
{
    boundaries_.controlAfterMove();
}



// control
void Foam::polyMoleculeCloud::controlBeforeForces()
{
    controllers_.controlBeforeForces();
}

void Foam::polyMoleculeCloud::clearLagrangianFields()
{
    iterator mol(this->begin());

    // Set accumulated quantities to zero
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().a() = vector::zero;

        mol().tau() = vector::zero;

        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;

        mol().R() = GREAT;
    }
}

void Foam::polyMoleculeCloud::calculateForce()
{
    calculatePairForces();
}

// update acceleration from net forces
void Foam::polyMoleculeCloud::updateAcceleration()
{
    accelerationUpdate();
}

void Foam::polyMoleculeCloud::accelerationUpdate()
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            mol().updateAcceleration(cP_);
        }
    }
}

// control
void Foam::polyMoleculeCloud::controlAfterForces()
{
    boundaries_.controlAfterForces();
    controllers_.controlAfterForces();
}

void Foam::polyMoleculeCloud::controlAfterVelocity()
{
    controllers_.controlVelocitiesII();    
}

void Foam::polyMoleculeCloud::postTimeStep()
{
    fields_.calculateFields();
    fields_.writeFields();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    controllers_.calculateStateProps();
    controllers_.outputStateResults();

    if(mesh_.time().outputTime())
    {
        writeReferredCloud();
    }

    trackingInfo_.clean(); 
}










//- used if you want to read a new field at every time-step from an input file
//- e.g. to be used in a utility that computes measurements
// Used by reconstructXmol utility - reconstructPar does not produce the XMOL
// files after parallel processing
void Foam::polyMoleculeCloud::readNewField()
{
    label initialSize = this->size();

    clear();

    IOPosition<Cloud<polyMolecule> > ioP(*this);

    if (ioP.headerOk())
    {
        ioP.readData(*this, false);
        ioP.close();
    }
    else
    {
        // WARNING
        WarningIn("readNewField()")
            << "Cannot read particle positions file " << nl
            << "    " << ioP.objectPath() << nl
            << "    assuming the initial cloud contains 0 particles." << endl;        
    }    
    
    particle::readFields(*this);

    polyMolecule::readFields(*this);

    if (this->size() != initialSize)
    {
        Info << "Changed polyMoleculeCloud size, from: " 
                << initialSize << ", to: " << this->size() << endl;
    }

    setSiteSizesAndPositions();    
}




void Foam::polyMoleculeCloud::setIPL()
{
   forAll(ipl_, c)
   {
        ipl_[c].clear();
   }

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        forAll(iL_.inverseDIL()[mol().cell()], c)
        {
            ipl_[iL_.inverseDIL()[mol().cell()][c]].append(&mol());
        }
    }
}

void Foam::polyMoleculeCloud::rebuildCellOccupancy()
{
    buildCellOccupancy();
}

void Foam::polyMoleculeCloud::prepareInteractions()
{
    iL_.setReferredParticles(cellOccupancy());
}

void Foam::polyMoleculeCloud::calculatePairForces()
{

    prepareInteractions();

    polyMolecule* molI = NULL;
    polyMolecule* molJ = NULL;

    {
        // Real-Real interactions
        const labelListList& dil = iL_.dil();

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

				forAll(dil[d], interactingCells)
				{
					List<polyMolecule*> cellJ =
						cellOccupancy_[dil[d][interactingCells]];

					forAll(cellJ, cellJMols)
					{
						molJ = cellJ[cellJMols];

						evaluatePair(molI, molJ);
					}
				}

				forAll(cellOccupancy_[d], cellIOtherMols)
				{
					molJ = cellOccupancy_[d][cellIOtherMols];

					if (molJ > molI)
					{
						evaluatePair(molI, molJ);
					}
				}
            }
        }
    }
    {
        // Real-Referred interactions
        forAll(iL_.refCellsParticles(), r)
        {
            const List<label>& realCells = iL_.refCells()[r].neighbouringCells();

            forAll(iL_.refCellsParticles()[r], i)
            {
            	molJ = iL_.refCellsParticles()[r][i];

				forAll(realCells, rC)
				{
					List<polyMolecule*> molsInCell = cellOccupancy_[realCells[rC]];

					forAll(molsInCell, j)
					{
						molI = molsInCell[j];
						evaluatePair(molI, molJ);
					}
				}
            }
        }
    }

}

void Foam::polyMoleculeCloud::writeXYZ(const fileName& fName) const
{
    OFstream os(fName);

    os << nSites() << nl << "polyMoleculeCloud site positions in angstroms" << nl;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
//         const polyMolecule::constantProperties& cP = constProps(mol().id());

        forAll(mol().sitePositions(), i)
        {
            const point& sP = mol().sitePositions()[i];

//             os << pot_.siteIdList()[cP.sites()[i].siteId()]
            os << cP_.siteNames(mol().id())[i]
                << ' ' << sP.x()*redUnits_.refLength()*1.0e10
                << ' ' << sP.y()*redUnits_.refLength()*1.0e10
                << ' ' << sP.z()*redUnits_.refLength()*1.0e10
                << nl;
        }
    }
}

// new function added to write the referred cloud 
// - to visualise the particles in ParaFOAM/VMD
void Foam::polyMoleculeCloud::writeReferredCloud()
{
    if(iL_.write())
    {
        const Time& runTime = mesh_.time();
        
        Info << "Writing out referred cloud" << endl;

        fileName timePath(runTime.path()/runTime.timeName()/"lagrangian");
        
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }
        
        fileName fName1(timePath/"referredCloud.xmol"); // VMD          
        fileName fName2(timePath/"referredCloud_RU.xmol");  //ParaFOAM
        
        label nParticles = iL_.referredCloud().size();
        
        Info << "number of particles = " << nParticles << endl;

        label nSites = 0;
        
        forAllIter
        (
            IDLList<polyMolecule>,
            iL_.referredCloud(),
            mol
        )
        {
            nSites += mol().sitePositions().size();
        }

        Info << "number of sites = " << nSites << endl;
        
        OFstream os1(fName1);
        OFstream os2(fName2);
        
        os1 << nSites << nl << "referred polyMoleculeCloud site positions in angstroms" << nl;
        os2 << nSites << nl << "referred polyMoleculeCloud site positions in reduced units" << nl;    
        
        forAllIter
        (
            IDLList<polyMolecule>,
            iL_.referredCloud(),
            mol
        )
        {
//             const polyMolecule::constantProperties& cP = constProps(mol().id());

            forAll(mol().sitePositions(), j)
            {            
            	const point& sP = mol().sitePositions()[j];

//                 os1 << pot_.siteIdList()[cP.sites()[j].siteId()]
                    os1 << cP_.siteNames(mol().id())[j]
                        << ' ' << sP.x()*redUnits_.refLength()*1e10
                        << ' ' << sP.y()*redUnits_.refLength()*1e10
                        << ' ' << sP.z()*redUnits_.refLength()*1e10
                        << nl;
                        
//                 os2 << pot_.siteIdList()[cP.sites()[j].siteId()]
                   os2 << cP_.siteNames(mol().id())[j]    
                        << ' ' << sP.x()
                        << ' ' << sP.y()
                        << ' ' << sP.z()
                        << nl;
            }
        }
    }
}

void Foam::polyMoleculeCloud::updateTrackingNumbersAfterRead()
{
    const_iterator mol(this->begin());
    
    label tN = 0;
    
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        if(mol().trackingNumber() > tN)
        {
            tN = mol().trackingNumber();
        }
    } 
    
    //- parallel-processing
    if(Pstream::parRun())
    {

        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << tN;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label trackingNumberProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> trackingNumberProc;
                }

                if(trackingNumberProc > tN)
                {
                    tN = trackingNumberProc;
                }
            }
        }
    }    
    
    moleculeTracking_.trackingIndex() = tN+1;
}

Foam::label Foam::polyMoleculeCloud::getTrackingNumber()
{
    return moleculeTracking_.getTrackingNumber();
}

// This function is not being used.It is a test to see if the max number of possible 
// labels have been exceeded. The maximum is usually so large that this is hardly ever possible.
// It may only be required for systems with an enormous turnover of molecules (adds and deletes).   
void Foam::polyMoleculeCloud::resetTrackingNumbers()
{
    moleculeTracking_.resetTrackingNumbers();

    if(moleculeTracking_.resetTracking())
    {
        iterator mol(this->begin());

        for (mol = this->begin(); mol != this->end(); ++mol)
        {
            mol().trackingNumber() = getTrackingNumber();
        }
    }
}



void Foam::polyMoleculeCloud::insertMolInCellOccupancy(polyMolecule* mol)
{
    cellOccupancy_[mol->cell()].append(mol);
}

void Foam::polyMoleculeCloud::removeMolFromCellOccupancy
(
    polyMolecule* molI
)
{
    DynamicList<polyMolecule*> updatedMolsInCell(0);

    const label& cellI = molI->cell();

    {
        const List<polyMolecule*>& molsInCell = cellOccupancy_[cellI];
    
        forAll(molsInCell, m)
        {
            polyMolecule* molJ = molsInCell[m];
    
            if(molI != molJ)
            {
                updatedMolsInCell.append(molJ);
            }
        }
    }

    cellOccupancy_[cellI].clear();
    cellOccupancy_[cellI].transfer(updatedMolsInCell);
}


void Foam::polyMoleculeCloud::removeMolFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<polyMolecule*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if(c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}



// ************************************************************************* //
