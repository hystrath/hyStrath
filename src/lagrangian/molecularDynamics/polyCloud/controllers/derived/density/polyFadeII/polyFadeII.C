/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    polyFadeII

Description

\*----------------------------------------------------------------------------*/

#include "polyFadeII.H"
#include "graph.H"

namespace Foam
{

void polyFadeII::readProperties()
{
    maxMolTries_ = 50;

    tauT_ = readScalar(propsDict_.lookup("tauT"));
    n_ = readLabel(propsDict_.lookup("polynomialDegree"));    
    
    if (propsDict_.found("empty"))
    {
        empty_ = Switch(propsDict_.lookup("empty"));
    }       
    
    velocity_ = vector::zero;
    temperature_ = 0.0;
    
    if (propsDict_.found("velocityOption"))
    {
        const word velocityOption = propsDict_.lookup("velocityOption");
        velocityOption_ = velocityOption;
        
        if(velocityOption_ == "maxwellian")
        {
            
        }
        else
        {
            FatalErrorIn("atomisticFadeII::atomisticFadeII()")
                << "Cannot find velocity option: " << velocityOption_
                << exit(FatalError);          
        }
    }   
    else
    {
        velocityOption_ = "maxwellian";
    }    
}


void polyFadeII::setVelocity
(
    polyMoleculeCloud& molCloud,
    vector& velocity
)
{
    if(velocityOption_ == "maxwellian")
    {
        velocity = velocity_ + getMaxwellianVelocity(molCloud);
    }
}

vector polyFadeII::getMaxwellianVelocity
(
    polyMoleculeCloud& molCloud
)
{
    const scalar& mass = molCloud.cP().mass(molId_);
    
    return sqrt(molCloud.redUnits().kB()*temperature_/mass)*vector
    (
        rndGen_.GaussNormalMD<scalar>(),
        rndGen_.GaussNormalMD<scalar>(),
        rndGen_.GaussNormalMD<scalar>()
    );
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
polyFadeII::polyFadeII()
:
    propsDict_(),
    rndGen_(clock::getTime(), 1),
    molId_(-1),
    nInserted_(0),
    nDeleted_(0),    
    boundedBox_(false),
    biggestRadius_(false),
    closestPositionBoundBox_(false),
    closestPositionCell_(false),
    exactPosition_(false),   
    empty_(false), 
    trackingNumbersIns_(0),
    trackingNumbersDel_(0)
{}

// Construct from cachedRandomMD
polyFadeII::polyFadeII
(
	const cachedRandomMD& rndMD
)
:
    propsDict_(),
// 	rndGen_(rndMD, false), //Create local copy of cached random values
    rndGen_(rndMD.seed(), rndMD.cacheSizeMult()),
    molId_(-1),
    nInserted_(0),
    nDeleted_(0),
    boundedBox_(false),
    biggestRadius_(false),
    closestPositionBoundBox_(false),
    closestPositionCell_(false),
    exactPosition_(false),
    empty_(false),
    trackingNumbersIns_(0),
    trackingNumbersDel_(0)
{}

// Construct from dict
polyFadeII::polyFadeII
(
    const dictionary& dict,
	const cachedRandomMD& rndMD
)
:
    propsDict_(dict.subDict("fadeProperties")),
// 	rndGen_(rndMD, false), //Create local copy of cached random values
    rndGen_(rndMD.seed(), rndMD.cacheSizeMult()),
	molId_(-1),
    nInserted_(0),
    nDeleted_(0),
    boundedBox_(false),
    biggestRadius_(false),
    closestPositionBoundBox_(false),
    closestPositionCell_(false),
    exactPosition_(false),   
    empty_(false), 
    trackingNumbersIns_(0),
    trackingNumbersDel_(0)
{
    timeIns_.clear();
    timeDel_.clear();
}

polyMolecule* polyFadeII::insertMoleculeInCloud
(
    polyMoleculeCloud& molCloud,
    const vector& position,
    const label& cell,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    const vector& specialPosition,
    const label molId,
    const label special
)
{
	label cellI = -1;
    label tetFaceI = -1;
    label tetPtI = -1;

    molCloud.mesh().findCellFacePt
    (
        position,
        cellI,
        tetFaceI,
        tetPtI
    );
                
    molCloud.createMolecule
    (
        position,
        cellI,
        tetFaceI,
        tetPtI,
        Q,
        v,
        a,
        pi,
        tau,
        specialPosition,
        special,
        molId,
        0.0,
        molCloud.getTrackingNumber()
    );

    polyMolecule* newMol = molCloud.last();

    molCloud.updateNeighbouringRadii(newMol);

    molCloud.insertMolInCellOccupancy(newMol);          
            
    return newMol;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyFadeII::~polyFadeII()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyFadeII::createInitialConfiguration
(
    const dictionary& dict,
    const label& molId,
    const scalar& deltaT,
    const word& insertOption,
    const word& deleteOption
)
{
    propsDict_ = dict.subDict("fadeProperties");
    molId_ = molId;
    
    deltaT_ = deltaT;
    
    readProperties();

    tauT_ -= deltaT_*5;
    
    trackingNumbersIns_.clear();
    trackingNumbersDel_.clear();
    timeIns_.clear();
    timeDel_.clear();    
    
    setInsertOption(insertOption);
    setDeleteOption(deleteOption);    
}

void polyFadeII::setInsertOption(const word& option)
{
    if(option == "biggestRadius")
    {
        biggestRadius_ = true;
        boundedBox_ = true;        
    }
/*    else if(option == "biggestRadiusII")
    {
        biggestRadiusII_ = true;
        boundedBox_ = true;        
    }   */ 
    else if(option == "closestPositionCell")
    {
        closestPositionCell_ = true;
    }
    else if(option == "closestPositionBoundBox")
    {
        closestPositionBoundBox_ = true;
        boundedBox_ = true;        
    }    
    else if(option == "exactPosition")
    {
        exactPosition_ = true;
    }
    else
    {
        FatalErrorIn("polyFadeII::polyFadeII()")
            << "Cannot find insertion option: " << option
            << exit(FatalError);    
    }
}

void polyFadeII::setDeleteOption(const word& option)
{
    deleteClosestPosition_ = false;
    deleteClosestPositionBoundBox_ = false;
    
    if(option == "closestPosition")
    {
        deleteClosestPosition_ = true;
    }
    else if(option == "closestPositionBoundBox")
    {
        deleteClosestPositionBoundBox_ = true;
    }    
    else
    {
        FatalErrorIn("polyFadeII::polyFadeII()")
            << "Cannot find deletion option: " << option
            << exit(FatalError);    
    }
}

void polyFadeII::controlMolecules
(
    polyMoleculeCloud& molCloud,
    const label& insertOrDelete,
    const boundedBox& bb,
    const List<vector>& positions 
)
{

    // reset properties
    nInserted_ = 0;
    nDeleted_ = 0;
    
    DynamicList<label> tNsIns(0);
    DynamicList<label> tNsDel(0);
    
    DynamicList<scalar> timeIns(0);
    DynamicList<scalar> timeDel(0);
    
    // find those molecules in the bound box
    
    DynamicList<polyMolecule*> molsInBB;
    
    if(boundedBox_)
    {
        IDLList<polyMolecule>::iterator mol(molCloud.begin());
        
        for
            (
             mol = molCloud.begin();
             mol != molCloud.end();
             ++mol
             )
        {
            if(bb.contains(mol().position()))
            {
                polyMolecule* molI = &mol();
                
                label tN = molI->trackingNumber();
                
                if(isNotFading(tN))
                {
                    molsInBB.append(molI);
                }
            }
        }
    }    
    
    //molsInBB.shrink();

    DynamicList<label> failedTNs(0);    
    
    // - insert in parallel

    for (int p = 0; p < Pstream::nProcs(); p++)
    {
        if(p == Pstream::myProcNo())
        {
            if(insertOrDelete > 0)
            {
                for(label i = 0; i < insertOrDelete; i++)
                {
                    label tN = -1;
                    
                    if(biggestRadius_)
                    {
                        tN = insertMoleculeInBoundBox
                        (
                            molCloud,
                            bb,
                            molsInBB,
                            failedTNs
                        );
                    }
                    else if(closestPositionBoundBox_)
                    {
                        tN = insertMoleculeAtClosestPositionBoundBox
                        (
                            molCloud,
                            bb,
                            molsInBB,
                            failedTNs,
                            positions[i]
                        );                    
                    }
                    else if(closestPositionCell_)
                    {
                        tN = insertMoleculeAtClosestPositionCell
                        (
                            molCloud,
                            bb,
                            failedTNs,                         
                            positions[i]
                        );                    
                    }                    
                    else if(exactPosition_)
                    {
                        tN = insertMoleculeAtExactPosition
                        (
                            molCloud,
                            bb,
                            positions[i]
                        );                    
                    }
                    
                    if(tN != -1)
                    {
                        tNsIns.append(tN);
                        timeIns.append(0.0);
                        nInserted_++;
                    }
                }
            }
        }

        molCloud.prepareInteractions();

        molCloud.updateRadii();
    }
    
    if(Pstream::parRun())
    {
        reduce(nInserted_, sumOp<label>()); 
    }    
    
    Info << "Inserted molecules (total): " << nInserted_ << endl;  

    // - delete
    
    if(insertOrDelete < 0)
    {
        for(label i = 0; i < mag(insertOrDelete); i++)
        {
            label tN = -1;
            
            if(deleteClosestPositionBoundBox_)
            {
                tN = deleteMoleculeFromClosestPositionBoundBox
                (
                    molCloud,
                    bb,
                    positions[i]
                );
            }
            else if(deleteClosestPosition_)  
            {
                tN = deleteMoleculeFromClosestPosition
                (
                    molCloud,
                    bb,
                    positions[i]
                );
            }            

            
//             Info << "Deleted tN = " << tN << endl;
            
            if(tN != -1)
            {
                nDeleted_++;
                timeDel.append(0.0);                
                tNsDel.append(tN);                
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(nDeleted_, sumOp<label>()); 
    }
    
    Info << "Deleted molecules (total): " << nDeleted_ << endl;

    insertInLists(tNsIns, tNsDel, timeIns, timeDel);
}




// void polyFadeII::deleteMoleculeFromCloud
// (
//     polyMoleculeCloud& molCloud,
//     polyMolecule& mol
// )
// {
// //     Info << "molecule potential energy (before) = " << mol.potentialEnergy() << endl; 
//     
//     molCloud.updateNeighbouringForces(&mol, true);
// 
//     DynamicList<polyMolecule*> refMols;    
//     DynamicList<label> refIds(0);
//     
//     molCloud.kernel().deleteParticle(&mol, refMols, refIds);    
// 
//     molCloud.updateNeighbouringForcesReferred(refMols, refIds, true); 
//     
// //     Info << "molecule potential energy (after) = " << mol.potentialEnergy() << endl;
// 
//     
//     
//     molCloud.removeMolFromCellOccupancy(&mol);
// 
//     
//     //- remove polyMolecule from cloud
//     molCloud.deleteParticle(mol);
// }


// greatest radius
// label polyFadeII::deleteMoleculeFromBoundBox
// (
//     polyMoleculeCloud& molCloud,
//     const boundedBox& bb
// )
// {
//    
//     label trackingNumber = -1;
//     scalar R = 0.0;
//     
//     {
//         IDLList<polyMolecule>::iterator mol(molCloud.begin());
// 
//         for
//         (
//             mol = molCloud.begin();
//             mol != molCloud.end();
//             ++mol
//         )
//         {
// //             if(bb.contains(mol().position()))
//             {
//                 polyMolecule* molI = &mol();
//                 
//                 if(molI->id() == molId_)
//                 {
//                     if(molI->R() > R)
//                     {
//                         R = molI->R();
//                         trackingNumber = molI->trackingNumber();
//                     }
//                 }
//             }
//         }
//     }
//     
//     return trackingNumber;
// }

label polyFadeII::deleteMoleculeFromClosestPositionBoundBox
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb, 
    const vector& r
)
{
   
    label trackingNumber = -1;
    scalar deltaR = GREAT;
//     vector position = vector::zero;
    
    {
        
        IDLList<polyMolecule>::iterator mol(molCloud.begin());

        for
        (
            mol = molCloud.begin();
            mol != molCloud.end();
            ++mol
        )
        {
            if(bb.contains(mol().position()))
            {
                polyMolecule* molI = &mol();
                
                if(molI->id() == molId_)
                {
                    scalar magRIJ = mag(molI->position() - r);
                    
                    if(magRIJ < deltaR)
                    {
                        deltaR = magRIJ;
                        trackingNumber = molI->trackingNumber();
//                         position = molI->position();
                    }
                }
            }
        }
    }
    
//     Info << "Deleting at position = " << position << endl;
    
    return trackingNumber;
}



label polyFadeII::deleteMoleculeFromClosestPosition
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb, 
    const vector& r
)
{
   
    label trackingNumber = -1;
    scalar deltaR = GREAT;
    
    {
        
        IDLList<polyMolecule>::iterator mol(molCloud.begin());

        for
        (
            mol = molCloud.begin();
            mol != molCloud.end();
            ++mol
        )
        {
            {
                polyMolecule* molI = &mol();
                
                if(molI->id() == molId_)
                {
                    scalar magRIJ = mag(molI->position() - r);
                    
                    if(magRIJ < deltaR)
                    {
                        deltaR = magRIJ;
                        trackingNumber = molI->trackingNumber();
                    }
                }
            }
        }
    }
    
    return trackingNumber;
}

void polyFadeII::deleteMolecules
(
    polyMoleculeCloud& molCloud,
    const DynamicList<label>& tNs
)
{
    DynamicList<polyMolecule*> molToDel;

    {
        IDLList<polyMolecule>::iterator mol(molCloud.begin());

        for
        (
            mol = molCloud.begin();
            mol != molCloud.end();
            ++mol
        )
        {
            if(findIndex(tNs, mol().trackingNumber()) != -1 )
            {
                polyMolecule* molI = &mol();
                molToDel.append(molI);
            }
        }
    }
    
    //molToDel.shrink();
    
    forAll(molToDel, i)
    {
        deleteMoleculeFromCloud(molCloud, *molToDel[i]);
    }
}

void polyFadeII::deleteMoleculeFromCloud
(
    polyMoleculeCloud& molCloud,
    polyMolecule& mol
)
{
    molCloud.removeMolFromCellOccupancy(&mol);

    //- remove polyMolecule from cloud
    molCloud.deleteParticle(mol);
}


void polyFadeII::checkFractions
(
    polyMoleculeCloud& molCloud
)
{
    DynamicList<label> tNsIns(0);
    DynamicList<label> tNsDel(0);
    
    label nInsDelMols = trackingNumbersIns_.size() + trackingNumbersDel_.size();    
    
    if(nInsDelMols > 0)
    {
        // update time first
        forAll(timeIns_, i)
        {
            timeIns_[i] += deltaT_;
        }
        forAll(timeDel_, i)
        {
            timeDel_[i] += deltaT_;
        }         
        
        IDLList<polyMolecule>::iterator mol(molCloud.begin());
    
        for
        (
            mol = molCloud.begin();
            mol != molCloud.end();
            ++mol
        )
        {
//             scalar frac = mol().fraction();

//             if(frac < 1.0)
            {
                label tN = mol().trackingNumber();
                
                // insertion
                label idIns = findIndex(trackingNumbersIns_, tN);
                
                if(idIns != -1)
                {
                    scalar t = timeIns_[idIns];

                    if(t < 0.5*tauT_)
                    {
                        mol().fraction() = 0.5*Foam::pow((2.0*t/tauT_), scalar(n_));
                    }
                    else
                    {
                        mol().fraction() = 1.0 
                                - 0.5*mag(Foam::pow((2.0*(t-tauT_)/tauT_), scalar(n_)));
                    }
                    
/*                    Info<< "time = " << t 
                        << ", insert fraction = " << mol().fraction()
                        << endl;*/
                    
                    if(mol().fraction() >= 1.0)
                    {
                        mol().fraction() = 1.0;
                        
                        tNsIns.append(tN);
                    }
                }
                
                // deletion
                label idDel = findIndex(trackingNumbersDel_, tN);
                
                if(idDel != -1)
                {
                    scalar t = timeDel_[idDel];

                    if(t < 0.5*tauT_)
                    {
                        mol().fraction() = 1.0 - 0.5*Foam::pow((2.0*t/tauT_), scalar(n_));
                    }
                    else
                    {
                        mol().fraction() = 
                                0.5*mag(Foam::pow((2.0*(t-tauT_)/tauT_), scalar(n_)));
                    }

//                     Info<< "time = " << t 
//                         << ", delete fraction = " << mol().fraction()
//                         << endl;

                    if(t >= tauT_)
                    {
                        mol().fraction() = 0;
                        tNsDel.append(tN);
                    }   
                }                  
            }
        }

        //tNsDel.shrink();
        //tNsIns.shrink();

        label sizeOfLists = tNsDel.size() + tNsIns.size();

        if(Pstream::parRun())
        {
            reduce(sizeOfLists, sumOp<label>());
        }

        if(sizeOfLists > 0)
        {
            deleteMolecules(molCloud, tNsDel);

            updateLists(tNsIns, tNsDel);
        }
        
        Info << "id: "<< molId_ << " no of mols inserting = " << trackingNumbersIns_.size()
            << ", no of mols deleting = " << trackingNumbersDel_.size()
            << endl;
    }
         
}



//- find a site and insert molecule close to another with the biggest radius R
label polyFadeII::insertMoleculeInBoundBox
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb,
    DynamicList<polyMolecule*>& molsInBB,
    DynamicList<label>& failedTNs
)
{    
    const polyMesh& mesh = molCloud.mesh();
    
    label trackingNumber = -1;
   
    //- choose an initial molecule - the one with the biggest radius R
    
    bool insertedMolecule = false;

    label nIter = 0;    
    
    while(!insertedMolecule && (nIter < maxMolTries_))
    {
        nIter++;
        
        label tN = -1;
        scalar rMax = 0.0;
        vector rI = vector::zero;     
        
        pickExistingMoleculeBiggestRadius(molsInBB, failedTNs, tN, rI, rMax);
        
        if(tN != -1)
        {
            vector rStart = vector::zero;
            label cellI = -1;
            
            if( chooseRandomPoint(mesh, bb, rI, rMax, cellI, rStart) )            
            {
                vector molVel(vector::zero);

                setVelocity(molCloud, molVel);

                scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                tensor Q = tensor
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                );   
                
                // insert molecule 
                polyMolecule* newMol = insertMoleculeInCloud
                (
                    molCloud,
                    rStart, //position
                    cellI,  //cell
                    Q, //Q
                    molVel, // velocity
                    vector::zero, // acc
                    vector::zero, // pi
                    vector::zero, // tau
                    vector::zero, // special position
                    molId_,
                    0 // special
                 );
                
                // tracking number
                trackingNumber = newMol->trackingNumber();  
                insertedMolecule = true;
            }

			failedTNs.append(tN);
			//failedTNs.shrink();
        }
        else if(empty_)
        {
            if(molsInBB.size() == 0)
            {
                Info << "Warning: empty bound box - trying to randomly insert" << endl;
                
                distributePoints box(bb, rndGen_); 
                
                vector rStart = box.randomPoint();

                label cellI = mesh.findCell(rStart);

                if(cellI != -1)
                {                    
                    vector molVel(vector::zero);
                    
                    setVelocity(molCloud, molVel);
                    
                    scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                    tensor Q = tensor
                    (
                        cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                        cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                        sin(psi)*sin(theta),
                        - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                        - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                        cos(psi)*sin(theta),
                        sin(theta)*sin(phi),
                        - sin(theta)*cos(phi),
                        cos(theta)
                    );   
                    
                    // insert molecule 
                    polyMolecule* newMol = insertMoleculeInCloud
                    (
                        molCloud,
                        rStart, //position
                        cellI,  //cell
                        Q, //Q
                        molVel, // velocity
                        vector::zero, // acc
                        vector::zero, // pi
                        vector::zero, // tau
                        vector::zero, // special position
                        molId_,
                        0 // special
                    );
                    
                    // tracking number
                    trackingNumber = newMol->trackingNumber();  
                    insertedMolecule = true;   
                }
            }
            else
            {
                nIter = maxMolTries_;
            }                
        }
        else
        {
            nIter = maxMolTries_;
        }
    }        

   
    return trackingNumber;
}

//- find a site and insert molecule close to a target reference point in a bound box 
// (optimal FADE, slightly costly)
label polyFadeII::insertMoleculeAtClosestPositionBoundBox
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb,
    DynamicList<polyMolecule*>& molsInBB,
    DynamicList<label>& failedTNs,
    const vector& r
)
{    
    const polyMesh& mesh = molCloud.mesh();
    
    label trackingNumber = -1;
    
    bool insertedMolecule = false;
    
    label nIter = 0;    
    
    while(!insertedMolecule && (nIter < maxMolTries_))
    {
        nIter++;
        
        label tN = -1;
        scalar R = 0.0;
        vector rI = vector::zero;
        
        pickExistingMoleculeClosestPositionBoundBox(molsInBB, failedTNs, tN, rI, r, R);
        
        if(tN != -1)
        {
            vector rStart = vector::zero;
            label cellI = -1;
            
            if( chooseRandomPoint(mesh, bb, rI, R, cellI, rStart) )            
            {
                vector molVel(vector::zero);

                setVelocity(molCloud, molVel);
                
                scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                tensor Q = tensor
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                ); 
                
                // insert molecule 
                polyMolecule* newMol = insertMoleculeInCloud
                (
                    molCloud,
                    rStart, //position
                    cellI,  //cell
                    Q, //Q
                    molVel, // velocity
                    vector::zero, // acc
                    vector::zero, // pi
                    vector::zero, // tau
                    vector::zero, // special position
                    molId_,
                    0 // special
                 );
                
                // tracking number
                trackingNumber = newMol->trackingNumber();  
                insertedMolecule = true;
            }

			failedTNs.append(tN);
			//failedTNs.shrink();
        }
        else if(empty_)
        {
            if(molsInBB.size() == 0)
            {
                Info << "Warning: empty bound box - inserting at defined position" << endl;
                
                vector rStart = r;
                
                label cellI = mesh.findCell(rStart);

                if(cellI != -1)
                {                    
                    vector molVel(vector::zero);
                    setVelocity(molCloud, molVel);
                    
                    scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                    tensor Q = tensor
                    (
                        cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                        cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                        sin(psi)*sin(theta),
                        - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                        - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                        cos(psi)*sin(theta),
                        sin(theta)*sin(phi),
                        - sin(theta)*cos(phi),
                        cos(theta)
                    ); 
                    
                    // insert molecule 
                    polyMolecule* newMol = insertMoleculeInCloud
                    (
                        molCloud,
                        rStart, //position
                        cellI,  //cell
                        Q, //Q
                        molVel, // velocity
                        vector::zero, // acc
                        vector::zero, // pi
                        vector::zero, // tau
                        vector::zero, // special position
                        molId_,
                        0 // special
                    );
                    
                    // tracking number
                    trackingNumber = newMol->trackingNumber();  
                    insertedMolecule = true;   
                }
            }
            else
            {
                nIter = maxMolTries_;
            }                
        }        
        
        else
        {
            nIter = maxMolTries_;
        }
    }        

   
    return trackingNumber;
}


// this picks a starting point from the closest molecule at r within its cell
// note this is still Optimal FADE, but it also accounts for empty cells
label polyFadeII::insertMoleculeAtClosestPositionCell
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb,
    DynamicList<label>& failedTNs,
    const vector& r
)
{    
    const polyMesh& mesh = molCloud.mesh();
    
    label trackingNumber = -1;
    
    bool insertedMolecule = false;
    
    label nIter = 0;    
    
    while(!insertedMolecule && (nIter < maxMolTries_))
    {
        nIter++;
        
        label tN = -1;
        scalar R = 0.0;
        vector rI = vector::zero;     

        pickExistingMoleculeClosestPositionCell(molCloud, failedTNs, tN, rI, r, R);
        
        if(tN != -1)
        {
            vector rStart = vector::zero;
            label cellI = -1;
            
            if( chooseRandomPoint(mesh, bb, rI, R, cellI, rStart) )            
            {
                vector molVel(vector::zero);
                setVelocity(molCloud, molVel);
                
                scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                tensor Q = tensor
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                ); 
                
                // insert molecule 
                polyMolecule* newMol = insertMoleculeInCloud
                (
                    molCloud,
                    rStart, //position
                    cellI,  //cell
                    Q, //Q
                    molVel, // velocity
                    vector::zero, // acc
                    vector::zero, // pi
                    vector::zero, // tau
                    vector::zero, // special position
                    molId_,
                    0 // special
                );
                
                // tracking number
                trackingNumber = newMol->trackingNumber();  
                insertedMolecule = true;
            }

            failedTNs.append(tN);
            //failedTNs.shrink();
        }
        else if(empty_)
        {
            label cellI = mesh.findCell(r);

            if(cellI != -1)
            {
                const List<polyMolecule*>& molsInCell = molCloud.cellOccupancy()[cellI];            
                
                if(molsInCell.size() == 0)
                {
                    Info << "Warning: empty cell - inserting at mid point" << endl;
                    
                    vector rStart = r;
                    
                    vector molVel(vector::zero);
                    setVelocity(molCloud, molVel);
                    
                    scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
                    scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

                    tensor Q = tensor
                    (
                        cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                        cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                        sin(psi)*sin(theta),
                        - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                        - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                        cos(psi)*sin(theta),
                        sin(theta)*sin(phi),
                        - sin(theta)*cos(phi),
                        cos(theta)
                    ); 
                    
                    // insert molecule 
                    polyMolecule* newMol = insertMoleculeInCloud
                    (
                        molCloud,
                        rStart, //position
                        cellI,  //cell
                        Q, //Q
                        molVel, // velocity
                        vector::zero, // acc
                        vector::zero, // pi
                        vector::zero, // tau
                        vector::zero, // special position
                        molId_,
                        0 // special
                    );
                    
                    // tracking number
                    trackingNumber = newMol->trackingNumber();  
                    insertedMolecule = true;   
                }
                else
                {
                    nIter = maxMolTries_;
                }                
            }
            else
            {
                nIter = maxMolTries_;
            }
        }
        else
        {
            nIter = maxMolTries_;
        }
    }
   
    return trackingNumber;
}


//- find a site and insert molecule close to a target reference point
label polyFadeII::insertMoleculeAtExactPosition
(
    polyMoleculeCloud& molCloud,
    const boundedBox& bb,
    const vector& r
)
{    
    const polyMesh& mesh = molCloud.mesh();
    
    label trackingNumber = -1;
    
    label cellI = mesh.findCell(r);
    
    if(cellI != -1)
    {
        vector molVel(vector::zero);
        setVelocity(molCloud, molVel);
        
        scalar phi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
        scalar theta(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);
        scalar psi(rndGen_.sample01<scalar>()*constant::mathematical::twoPi);

        tensor Q = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        ); 
        
        // insert molecule 
        polyMolecule* newMol = insertMoleculeInCloud
        (
            molCloud,
            r, //position
            cellI,  //cell
            Q, //Q
            molVel, // velocity
            vector::zero, // acc
            vector::zero, // pi
            vector::zero, // tau
            vector::zero, // special position
            molId_,
            0 // special
        );
        
        // tracking number
        trackingNumber = newMol->trackingNumber();  
    }
    else
    {
        Info << "ERROR - > TRYING TO INSERT MOLECULE OUTSIDE MESH " << endl;   
    }
    
    return trackingNumber;
}


void polyFadeII::pickExistingMoleculeBiggestRadius
(
     const DynamicList<polyMolecule*>& molsInBB,
     const DynamicList<label>& failedTNs,
     label& tN,
     vector& rI,
     scalar& R
)
{

    forAll(molsInBB, i)
    {
        polyMolecule* molI = molsInBB[i];

		if(findIndex(failedTNs, molI->trackingNumber()) == -1)
		{
			if(molI->R() > R)
			{
				R = molI->R();
				tN = molI->trackingNumber();
				rI = molI->position();
			}
		}
    }
}

// this function is more computational efficient, but not general 
// (for e.g. the FENE case, requires another option)
void polyFadeII::pickExistingMoleculeClosestPositionCell
(
    polyMoleculeCloud& molCloud,
    const DynamicList<label>& failedTNs,
    label& tN,
    vector& rI,
    const vector& r,
    scalar& R     
)  
{
   const polyMesh& mesh = molCloud.mesh();
    
    scalar deltaR = GREAT;
    
    label cellI = mesh.findCell(r);
    
    if(cellI != -1)
    {
        const List<polyMolecule*>& molsInCell = molCloud.cellOccupancy()[cellI];

        forAll(molsInCell, i)
        {
            polyMolecule* molI = molsInCell[i];

			if(findIndex(failedTNs, molI->trackingNumber()) == -1)
			{
				scalar magRIJ = mag(molI->position() - r);

				if(magRIJ < deltaR)
				{
					deltaR = magRIJ;
					tN = molI->trackingNumber();
					rI = molI->position();
					R = molI->R();
				}
			}
        } 
    }    
}

void polyFadeII::pickExistingMoleculeClosestPositionBoundBox
(
    const DynamicList<polyMolecule*>& molsInBB,
    const DynamicList<label>& failedTNs,
    label& tN,
    vector& rI,
    const vector& r,
    scalar& R     
)  
{
    
    scalar deltaR = GREAT;

    forAll(molsInBB, i)
    {
        polyMolecule* molI = molsInBB[i];
        
		if(findIndex(failedTNs, molI->trackingNumber()) == -1)
		{
			scalar magRIJ = mag(molI->position() - r);

			if(magRIJ < deltaR)
			{
				deltaR = magRIJ;
				tN = molI->trackingNumber();
				rI = molI->position();
				R = molI->R();
			}
		}
    }    
}


bool polyFadeII::chooseRandomPoint
(
    const polyMesh& mesh,
    const boundedBox& bb,
    const vector& rI, 
    const scalar& R,
    label& cellI,
    vector& rStart
)
{
    bool chosenPoint = false;
    
    label nIter = 0;
    
    while ( !chosenPoint && (nIter < 40) )
    {   
        nIter++;
        
        //- select a random direction
        scalar magV = 0.0;
        vector randDirection(vector::zero);

        while ( !(magV > 0.0) )
        {
        	randDirection = rndGen_.sampleVectorMD<vector>();
        	magV = mag(randDirection);
        }

        //- normalise the random vector (unit vector)
        randDirection /= mag(randDirection);

        //- the starting point is selected to be the point in which it is most
        //  likely to find a gap

        rStart = rI + (0.5*R*randDirection);
        
        if(bb.contains(rStart))
        {
            cellI = mesh.findCell(rStart);
            
            if(cellI != -1)
            {
                chosenPoint = true;
            }
        }
    }
    
    return chosenPoint;
}


void polyFadeII::updateLists
(
    DynamicList<label>& tNsIns,
    DynamicList<label>& tNsDel
)
{
    //tNsIns.shrink();
    //tNsDel.shrink();

    // sync processors 
    if (Pstream::parRun())
    {
        List<label> tNsInsList(tNsIns.size(), 0);
        List<label> tNsDelList(tNsDel.size(), 0);

        forAll(tNsIns, i)
        {
            tNsInsList[i]=tNsIns[i];
        }

        forAll(tNsDel, i)
        {
            tNsDelList[i]=tNsDel[i];
        }

        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << tNsInsList << tNsDelList;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<label> tNsInsListProc;
                List<label> tNsDelListProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> tNsInsListProc >> tNsDelListProc;
                }
                
                forAll(tNsInsListProc, i)
                {
                    tNsIns.append(tNsInsListProc[i]);
                }

                forAll(tNsDelListProc, i)
                {
                    tNsDel.append(tNsDelListProc[i]);
                }             
            }
        }
    }
    
    DynamicList<label> tNsInsNew(0);
    DynamicList<label> tNsDelNew(0);    
    
    DynamicList<scalar> timeInsNew(0);
    DynamicList<scalar> timeDelNew(0);
    
    forAll(trackingNumbersIns_, i)
    {
        if(findIndex(tNsIns, trackingNumbersIns_[i]) == -1)
        {
            tNsInsNew.append(trackingNumbersIns_[i]);
            timeInsNew.append(timeIns_[i]);
        }
    }

    forAll(trackingNumbersDel_, i)
    {
        if(findIndex(tNsDel, trackingNumbersDel_[i]) == -1)
        {
            tNsDelNew.append(trackingNumbersDel_[i]);
            timeDelNew.append(timeDel_[i]);
        }
    }    

    trackingNumbersIns_.clear();
    trackingNumbersDel_.clear();
    
    timeIns_.clear();
    timeDel_.clear();
    
    //trackingNumbersIns_.transfer(tNsInsNew.shrink());
    //trackingNumbersDel_.transfer(tNsDelNew.shrink());

    trackingNumbersIns_.transfer(tNsInsNew);
    trackingNumbersDel_.transfer(tNsDelNew);
    
    //timeIns_.transfer(timeInsNew.shrink());
    //timeDel_.transfer(timeDelNew.shrink());

    timeIns_.transfer(timeInsNew);
    timeDel_.transfer(timeDelNew);
}

void polyFadeII::insertInLists
(
    DynamicList<label>& tNsIns,
    DynamicList<label>& tNsDel,
    DynamicList<scalar>& timeIns,
    DynamicList<scalar>& timeDel     
)
{
    //tNsIns.shrink();
    //tNsDel.shrink();
    
    //timeIns.shrink();
    //timeDel.shrink();
    
    // sync processors 
    if (Pstream::parRun())
    {
        List<label> tNsInsList(tNsIns.size(), 0);
        List<label> tNsDelList(tNsDel.size(), 0);

        List<scalar> timeInsList(timeIns.size(), 0.0);
        List<scalar> timeDelList(timeDel.size(), 0.0);        
        
        forAll(tNsIns, i)
        {
            tNsInsList[i]=tNsIns[i];
            timeInsList[i]=timeIns[i];
        }

        forAll(tNsDel, i)
        {
            tNsDelList[i]=tNsDel[i];
            timeDelList[i]=timeDel[i];
        }

        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << tNsInsList << tNsDelList 
                                << timeInsList << timeDelList;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<label> tNsInsListProc;
                List<label> tNsDelListProc;
                List<scalar> timeInsListProc;
                List<scalar> timeDelListProc;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> tNsInsListProc >> tNsDelListProc
                                >> timeInsListProc >> timeDelListProc;
                }

                forAll(tNsInsListProc, i)
                {
                    trackingNumbersIns_.append(tNsInsListProc[i]);
                    timeIns_.append(timeInsListProc[i]);
                }
                

                forAll(tNsDelListProc, i)
                {
                    trackingNumbersDel_.append(tNsDelListProc[i]);
                    timeDel_.append(timeDelListProc[i]);
                }                
            }
        }
    }
    
    
    forAll(tNsIns, i)
    {
        trackingNumbersIns_.append(tNsIns[i]);
        timeIns_.append(timeIns[i]);
    }

    forAll(tNsDel, i)
    {
        trackingNumbersDel_.append(tNsDel[i]);
        timeDel_.append(timeDel[i]);
    }     
    
    //timeIns_.shrink();
    //timeDel_.shrink();
    //trackingNumbersIns_.shrink();
    //trackingNumbersDel_.shrink();
}

// update properties 
void polyFadeII::updateProperties
(
    const dictionary& newDict
)
{}

void polyFadeII::write
(
    const Time& time,
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyFadeII::output(Time& time)
{
    label nBins = label((tauT_/deltaT_)+0.5);
    
    scalarField timeField(nBins, 0.0);    
    scalarField fractionIns(nBins, 0.0);
    scalarField fractionDel(nBins, 0.0);    
    
    forAll(timeField, i)
    {
        timeField[i] = i*deltaT_;
        scalar t = timeField[i];

        if(t < 0.5*tauT_)
        {
            fractionIns[i] = 0.5*Foam::pow((2.0*t/tauT_), scalar(n_));
        }
        else
        {
            fractionIns[i] = 1.0 
                    - 0.5*mag(Foam::pow((2.0*(t-tauT_)/tauT_), scalar(n_)));
        }

        if(t < 0.5*tauT_)
        {
            fractionDel[i] = 1.0 - 0.5*Foam::pow((2.0*t/tauT_), scalar(n_));
        }
        else
        {
            fractionDel[i] = 
                    0.5*mag(Foam::pow((2.0*(t-tauT_)/tauT_), scalar(n_)));
        }
                    
        if(t >= tauT_)
        {
            fractionDel[i] = 0;
        }   
    }
    
    fileName casePath(time.path());
    
    writeTimeData
    (
        casePath,
        "fade_insert.xy",
        timeField,
        fractionIns
    ); 
    
    writeTimeData
    (
        casePath,
        "fade_delete.xy",
        timeField,
        fractionDel
    );     
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
