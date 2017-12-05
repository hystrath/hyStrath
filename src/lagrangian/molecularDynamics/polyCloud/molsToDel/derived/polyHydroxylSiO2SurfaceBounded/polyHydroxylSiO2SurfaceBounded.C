/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "polyHydroxylSiO2SurfaceBounded.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyHydroxylSiO2SurfaceBounded, 0);

addToRunTimeSelectionTable(polyMolsToDeleteModel, polyHydroxylSiO2SurfaceBounded, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void polyHydroxylSiO2SurfaceBounded::readInCylinder()
{
    dictionary propsDict(propsDict_.subDict("punchOutCylinderProperties"));
    radius_ = readScalar(propsDict.lookup("radius"));
    startPoint_ = propsDict.lookup("startPoint");
    endPoint_ = propsDict.lookup("endPoint");
    unitVector_ = (endPoint_ - startPoint_)/mag(endPoint_ - startPoint_);
    deltaR_ = readScalar(propsDict.lookup("deltaR")); 
    
}

void polyHydroxylSiO2SurfaceBounded::readInBoundBox()
{
    dictionary propsDict(propsDict_.subDict("boundBoxTopSurface"));
    
    vector startPoint = propsDict.lookup("startPoint");
    vector endPoint = propsDict.lookup("endPoint");
    hDirection_ = propsDict.lookup("hDirection");
    hDirection_ /= molCloud_.redUnits().refLength();
    checkBoundBox(box_, startPoint, endPoint);
}

// Construct from components
polyHydroxylSiO2SurfaceBounded::polyHydroxylSiO2SurfaceBounded
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMolsToDeleteModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    b_(readScalar(propsDict_.lookup("bondRadius"))),
   
    temperature_(readScalar(propsDict_.lookup("temperature"))),
    bulkVelocity_(propsDict_.lookup("bulkVelocity"))
{
    readInCylinder();
    readInBoundBox();
    
  
    
    // bond radius
    b_ /= molCloud_.redUnits().refLength();
        
    Info << " max bond radius = " << b_ << endl;    

    // set molecule ids

    const word siName = propsDict_.lookup("siliconIdName");
    const word oName = propsDict_.lookup("oxygenIdName");
    const word hName = propsDict_.lookup("hydrogenIdName");    

    const List<word>& idList(molCloud_.pot().idList());
    
   
    forAll(idList, i)
    {
        if(siName == idList[i])
        {
            SiId_ = i;
        }
        
        if(oName == idList[i])
        {
            OId_ = i;
        }       
        
        if(hName == idList[i])
        {
            HId_ = i;
        } 
    }
    
    Info<< " silicon id = " << SiId_ 
        << ", oxygen id = " << OId_
        << ", hydrogen id = " << HId_
        << endl;

    
    punchOutMolecules();
    
    check();
    
    checkForMoleculesInBox();
    
    deleteMoleculesInBox();
    
    check();
    
    checkForMoleculesInCylinder();
    
    deleteMoleculesInCylinder();
    
    check();
    
    checkForMoleculesInBox();  
    
//     addMoleculesInBox();
    
    addMoleculesInBoxAndCylinder();
    
    check();
    
    checkForMoleculesInBox();
    
    checkForMoleculesInCylinder();    
}


void polyHydroxylSiO2SurfaceBounded::punchOutMolecules()
{
    Info << nl << "Punching out cylinder " << endl; 
    
    DynamicList<polyMolecule*> molsToDel;
    
    label initialSize = molCloud_.size();    
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        scalar rSEMag = mag(endPoint_ - startPoint_);

        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            const vector& rI = mol().position();
            vector rSI = rI - startPoint_;
            scalar centreLineDistance = rSI & unitVector_;

            //- step 1: test polyMolecule is between starting point and end point
            if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
            {
                vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

                //step 2: test polyMolecule is within radial distance of centre-line
                if(mag(pointOnCentreLine-rI) <= radius_)
                {
                    label molId = mol().id();
                    
                    if( (molId == SiId_) || (molId == OId_) || (molId == HId_) )
                    {
                        polyMolecule* molI = &mol();
                        molsToDel.append(molI);
                    }
                }
            }
        }
    }
    
   // molsToDel.shrink();

    forAll(molsToDel, m)
    {
        deleteMolFromMoleculeCloud(*molsToDel[m]);
    }

    label molsKept = initialSize - molsToDel.size();

    Info<< tab << " initial polyMolecules: " <<  initialSize 
        << ", polyMolecules kept: " <<  molsKept
        << ", polyMolecules removed: " << molsToDel.size() 
        << endl;


    // as a precaution: rebuild cell occupancy
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareKernel();
}

void polyHydroxylSiO2SurfaceBounded::check()
{
    Info << nl << "Separating and counting atoms " << endl; 

    siliconAtoms_.clear();
    oxygenAtoms_.clear();
 
    {
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            if(molI().id() == SiId_)
            {
                siliconAtoms_.append(molI().trackingNumber());
            }
            else if(molI().id() == OId_)
            {
                oxygenAtoms_.append(molI().trackingNumber());
            }            
        }
    }
    
    //siliconAtoms_.shrink();
    //oxygenAtoms_.shrink();
    
    Info << " number of silicon atoms = " << siliconAtoms_.size() << endl;
    Info << " number of oxygen atoms = " << oxygenAtoms_.size() << endl;

    
    siliconToOxygenBonds_.clear();
    oxygenToSiliconBonds_.clear();
    siliconToOxygenBonds_.setSize(siliconAtoms_.size(), 0);
    oxygenToSiliconBonds_.setSize(oxygenAtoms_.size(), 0);   
    
    Info << nl << " checking for bonds" << endl; 
    
    
    {    
        label nMolsInt = 0;
        label nMolsExt = 0;
        label nMolsRef = 0;

        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        {
            polyMolecule* molI = &mol();
            polyMolecule* molJ = &mol();

            forAll(molCloud_.il(), c)
            {
                nMolsExt = molCloud_.il()[c].size();
                nMolsInt = molCloud_.cellOccupancy()[c].size();
                nMolsRef = molCloud_.kernel().referredInteractionList()[c].size();

                for (int i = 0; i < nMolsInt; i++)
                {
                    molI = molCloud_.cellOccupancy()[c][i];

                    for (int j = 0; j < nMolsInt; j++)
                    {
                        if(j > i)
                        {
                            molJ = molCloud_.cellOccupancy()[c][j];
                            
                            label siTN = molI->trackingNumber();
                            label oTN = molJ->trackingNumber();
                            
                            if( testForPair(molI->id(), molJ->id(), siTN, oTN) )
                            {
                                scalar rIJmag = mag(molI->position() - molJ->position());
                                
                                if(rIJmag <= b_)
                                {
                                    label ids = findIndex(siliconAtoms_, siTN);
                                    
                                    if(ids != -1)
                                    {
                                        siliconToOxygenBonds_[ids]++;
                                    }
                                    label ido = findIndex(oxygenAtoms_, oTN);
                                    
                                    if(ido != -1)
                                    {
                                        oxygenToSiliconBonds_[ido]++;
                                    }
                                }
                            }
                        }
                    }

                    for (int j = 0; j < nMolsExt; j++)
                    {
                        molJ = molCloud_.il()[c][j];
                        
                        label siTN = molI->trackingNumber();
                        label oTN = molJ->trackingNumber();
                        
                        if( testForPair(molI->id(), molJ->id(), siTN, oTN) )
                        {
                            scalar rIJmag = mag(molI->position() - molJ->position());
                            
                            if(rIJmag <= b_)
                            {
                                label ids = findIndex(siliconAtoms_, siTN);
                                
                                if(ids != -1)
                                {
                                    siliconToOxygenBonds_[ids]++;
                                }
                                label ido = findIndex(oxygenAtoms_, oTN);
                                
                                if(ido != -1)
                                {
                                    oxygenToSiliconBonds_[ido]++;
                                }
                            }
                        }
                    }

                    for (int j = 0; j < nMolsRef; j++)
                    {
                        molJ = molCloud_.kernel().referredInteractionList()[c][j];
                        
                        label siTN = molI->trackingNumber();
                        label oTN = molJ->trackingNumber();
                        
                        if( testForPair(molI->id(), molJ->id(), siTN, oTN) )
                        {
                            scalar rIJmag = mag(molI->position() - molJ->position());
                            
                            if(rIJmag <= b_)
                            {
                                label ids = findIndex(siliconAtoms_, siTN);
                                
                                if(ids != -1)
                                {
                                    siliconToOxygenBonds_[ids]++;
                                }
                                label ido = findIndex(oxygenAtoms_, oTN);
                                
                                if(ido != -1)
                                {
                                    oxygenToSiliconBonds_[ido]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // oxygen
    
    List<label> oBonds(10,0);
    List<label> sBonds(10,0);    
    
    forAll(oxygenToSiliconBonds_, i)
    {
        const label& nOB = oxygenToSiliconBonds_[i];
        oBonds[nOB]++;
    }
    
    forAll(siliconToOxygenBonds_, i)
    {
        const label& nsB = siliconToOxygenBonds_[i];
        sBonds[nsB]++;
    }
/*    Info << " system bond info: " << endl;    
    Info << " oxygenToSiliconBonds_ " << oBonds << endl;
    Info << " siliconToOxygenBonds_ " << sBonds << endl; */   

}

void polyHydroxylSiO2SurfaceBounded::checkForMoleculesInBox()
{
    // bound box details
    {
        List<label> oBondsBB(10,0);
        List<label> sBondsBB(10,0); 
        
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            if(box_.contains(molI().position()))
            {
                if(molI().id() == SiId_)
                {
                    label id = findIndex(siliconAtoms_,molI().trackingNumber());
                    const label& nsB = siliconToOxygenBonds_[id];
                    sBondsBB[nsB]++;
                }
                else if(molI().id() == OId_)
                {
                    label id = findIndex(oxygenAtoms_,molI().trackingNumber());
                    const label& noB = oxygenToSiliconBonds_[id];
                    oBondsBB[noB]++;
                }            
            }
        }
        
        Info << nl << "Information in bound-box "<< endl;

        Info << " oxygenToSiliconBonds_ " << oBondsBB << endl;
        Info << " siliconToOxygenBonds_ " << sBondsBB << endl;        
    }
}

void polyHydroxylSiO2SurfaceBounded::checkForMoleculesInCylinder()
{
    // bound box details
    {
        List<label> oBondsBB(10,0);
        List<label> sBondsBB(10,0); 
        
       IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        scalar rSEMag = mag(endPoint_ - startPoint_);

        for
        (
            mol = molCloud_.begin();
            mol != molCloud_.end();
            ++mol
        )
        {
            const vector& rI = mol().position();
            vector rSI = rI - startPoint_;
            scalar centreLineDistance = rSI & unitVector_;

            //- step 1: test polyMolecule is between starting point and end point
            if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
            {
                vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

                //step 2: test polyMolecule is within radial distance of centre-line
                if(mag(pointOnCentreLine-rI) <= (radius_ + deltaR_))
                {

                    if(mol().id() == SiId_)
                    {
                        label id = findIndex(siliconAtoms_, mol().trackingNumber());
                        const label& nsB = siliconToOxygenBonds_[id];
                        sBondsBB[nsB]++;
                    }
                    else if(mol().id() == OId_)
                    {
                        label id = findIndex(oxygenAtoms_,mol().trackingNumber());
                        const label& noB = oxygenToSiliconBonds_[id];
                        oBondsBB[noB]++;
                    }            
                }
            }
        }
        
        Info << nl << "Information in bound-box "<< endl;

        Info << " oxygenToSiliconBonds_ " << oBondsBB << endl;
        Info << " siliconToOxygenBonds_ " << sBondsBB << endl;        
    }
}

void polyHydroxylSiO2SurfaceBounded::deleteMoleculesInBox()
{
    DynamicList<polyMolecule*> molsToDel;
    
    {
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            if(box_.contains(molI().position()))
            {
                if(molI().id() == SiId_)
                {
                    label id = findIndex(siliconAtoms_, molI().trackingNumber());
                    
                    if(siliconToOxygenBonds_[id] < 4)
                    {
                        polyMolecule* mol = &molI();
                        molsToDel.append(mol);                        
                    }
                }
            }
        }
    }
        
    //molsToDel.shrink();
    
    Info << nl<< "Silicon molecules to delete = " << molsToDel.size() << endl;
    
    
    forAll(molsToDel, m)
    {
//         Info << "position = " << molsToDel[m]->position() << endl;
        
        molCloud_.deleteParticle(*molsToDel[m]);
    }    
    
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareKernel();
}

void polyHydroxylSiO2SurfaceBounded::deleteMoleculesInCylinder()
{
    DynamicList<polyMolecule*> molsToDel;
    
    {
       
       IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        scalar rSEMag = mag(endPoint_ - startPoint_);

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            const vector& rI = molI().position();
            vector rSI = rI - startPoint_;
            scalar centreLineDistance = rSI & unitVector_;

            //- step 1: test polyMolecule is between starting point and end point
            if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
            {
                vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

                //step 2: test polyMolecule is within radial distance of centre-line
                if(mag(pointOnCentreLine-rI) <= (radius_ + deltaR_))
                {
                    if(molI().id() == SiId_)
                    {
                        label id = findIndex(siliconAtoms_, molI().trackingNumber());
                        
                        if(siliconToOxygenBonds_[id] < 4)
                        {
                            polyMolecule* mol = &molI();
                            molsToDel.append(mol);                        
                        }
                    }
                }
            }
        }
    }
        
    //molsToDel.shrink();
    
    Info << nl<< "Silicon molecules to delete = " << molsToDel.size() << endl;
    
    
    forAll(molsToDel, m)
    {
//         Info << "position = " << molsToDel[m]->position() << endl;
        
        molCloud_.deleteParticle(*molsToDel[m]);
    }    
    
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareKernel();
}


void polyHydroxylSiO2SurfaceBounded::addMoleculesInBoxAndCylinder()
{
    DynamicList<polyMolecule*> molsToAddBox;
    DynamicList<polyMolecule*> molsToAddCylinder;
  
    {
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        scalar rSEMag = mag(endPoint_ - startPoint_);
        
        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {

            const vector& rI = molI().position();
            vector rSI = rI - startPoint_;
            scalar centreLineDistance = rSI & unitVector_;
            
            bool added = false;
            
            //- step 1: test polyMolecule is between starting point and end point
            if((centreLineDistance <= rSEMag) && (centreLineDistance >= 0.0))
            {
                vector pointOnCentreLine = centreLineDistance*unitVector_ + startPoint_;

                //step 2: test polyMolecule is within radial distance of centre-line
                if(mag(pointOnCentreLine-rI) <= (radius_ + deltaR_))
                {      
                    if(molI().id() == OId_)
                    {
                        label id = findIndex(oxygenAtoms_, molI().trackingNumber());
                        
                        if(oxygenToSiliconBonds_[id] < 2)
                        {
                            polyMolecule* mol = &molI();
                            molsToAddCylinder.append(mol); 
                            added = true;
                        }
                    }   
                }
            }
            
            if(!added)
            {    
                if(box_.contains(molI().position()))
                {
                    if(molI().id() == OId_)
                    {
                        label id = findIndex(oxygenAtoms_, molI().trackingNumber());
                        
                        if(oxygenToSiliconBonds_[id] < 2)
                        {
                            polyMolecule* mol = &molI();
                            molsToAddBox.append(mol); 
                        }
                    }            
                }
            }
        }
    }
    
    //molsToAddBox.shrink();
    //molsToAddCylinder.shrink();
    
    // add hydrogen molecules 
    Info << nl<< "Hydrogen molecules to add in box = " << molsToAddBox.size() << endl;    
    
    forAll(molsToAddBox, i)
    {
        vector newPosition = molsToAddBox[i]->position() + hDirection_;
        
//         Info << "position = " << newPosition << endl;
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            newPosition,
            cell,
            tetFace,
            tetPt
        );

        if(cell != -1)
        {
            insertMolecule
            (
                newPosition,
                cell,
                tetFace,
                tetPt,
                HId_,
                false,
                true,
                temperature_,
                bulkVelocity_
            );
        }
    }
    
    
    Info << nl<< "Hydrogen molecules to add in cylinder = " << molsToAddCylinder.size() << endl;
    
    forAll(molsToAddCylinder, i)
    {
        const vector&  rI = molsToAddCylinder[i]->position();
        
        vector rICentre = (
                            (((rI - startPoint_ ) & unitVector_)*unitVector_)
                            + startPoint_
                          ) - rI;
        
        rICentre /= mag(rICentre);
        
        vector newPosition = rI + mag(hDirection_)*rICentre;
        
//         Info << "position = " << newPosition << endl;
        
        label cell = mesh_.findCell(newPosition);

        if(cell != -1)
        {
            insertMolecule
            (
                newPosition,
                cell,
                HId_,
                false,
                true,
                temperature_,
                bulkVelocity_
            );
        }
    }    
    
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareKernel();    
}
/*
void polyHydroxylSiO2SurfaceBounded::addMoleculesInBox()
{
    DynamicList<polyMolecule*> molsToAdd;
    
    {
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            if(box_.contains(molI().position()))
            {
                if(molI().id() == OId_)
                {
                    label id = findIndex(oxygenAtoms_, molI().trackingNumber());
                    
                    if(oxygenToSiliconBonds_[id] < 2)
                    {
                        polyMolecule* mol = &molI();
                        molsToAdd.append(mol);                        
                    }
                }            
            }
        }
    }
    
    molsToAdd.shrink();
    
    // add hydrogen molecules 
    Info << nl<< "Hydrogen molecules to add = " << molsToAdd.size() << endl;    
    
    forAll(molsToAdd, i)
    {
        vector newPosition = molsToAdd[i]->position() + hDirection_;
        
//         Info << "position = " << newPosition << endl;
        
        label cell = mesh_.findCell(newPosition);

        if(cell != -1)
        {
            insertMolecule
            (
                newPosition,
                cell,
                HId_,
                false,
                true,
                temperature_,
                bulkVelocity_
            );
        }
    }
    
    molCloud_.rebuildCellOccupancy();
    molCloud_.prepareKernel();    
}*/

bool polyHydroxylSiO2SurfaceBounded::testForPair
(
    const label& idI,
    const label& idJ,
    label& tnI,
    label& tnJ
)
{
    bool correctPair = false;
    
    if(idI != idJ)
    {
        if((idI == SiId_) && (idJ == OId_))
        {
            correctPair = true;
        }
        
        if((idJ == SiId_) && (idI == OId_))
        {
            label tnTemp = tnI;
            tnI = tnJ;
            tnJ = tnTemp;
            
            correctPair = true;
        }        
        
    }
    
    return correctPair;
}

void polyHydroxylSiO2SurfaceBounded::checkBoundBox
(
    boundBox& b,
    const vector& startPoint,
    const vector& endPoint
)
{
    vector& vMin = b.min();
    vector& vMax = b.max();

    if(startPoint.x() < endPoint.x())
    {
        vMin.x() = startPoint.x();
        vMax.x() = endPoint.x();
    }
    else
    {
        vMin.x() = endPoint.x();
        vMax.x() = startPoint.x();
    }
    if(startPoint.y() < endPoint.y())
    {
        vMin.y() = startPoint.y();
        vMax.y() = endPoint.y();
    }
    else
    {
        vMin.y() = endPoint.y();
        vMax.y() = startPoint.y();
    }
    if(startPoint.z() < endPoint.z())
    {
        vMin.z() = startPoint.z();
        vMax.z() = endPoint.z();
    }
    else
    {
        vMin.z() = endPoint.z();
        vMax.z() = startPoint.z();
    }
}
/*
void polyHydroxylSiO2SurfaceBounded::insertMolecule
(
    const point& position,
    label& cell,
    const label& id,
    const bool& tethered,
    const bool& frozen,
    const scalar& temperature,
    const vector& bulkVelocity
)
{
    if (cell == -1)
    {
        cell = mesh_.findCell(position);
    }

    if (cell == -1)
    {
        FatalErrorIn("Foam::polyLatticeZone::insertMolecule()")
            << "Position specified does not correspond to a mesh cell." << nl
            << abort(FatalError);
    }

    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_TETHERED;
    }

    if (frozen)
    {
        specialPosition = position;

        special = polyMolecule::SPECIAL_FROZEN;
    }

    const polyMolecule::constantProperties& cP = molCloud_.constProps(id);
    
//     Info << "mass = " << cP.mass() << endl;
    
    vector v = equipartitionLinearVelocity(temperature, cP.mass());

    v += bulkVelocity;

    vector pi = vector::zero;

    tensor Q = I;

    if (!cP.pointMolecule())
    {
        pi = equipartitionAngularMomentum(temperature, id);

        scalar phi(molCloud_.rndGen().scalar01()*constant::mathematical::twoPi);

        scalar theta(molCloud_.rndGen().scalar01()*constant::mathematical::twoPi);

        scalar psi(molCloud_.rndGen().scalar01()*constant::mathematical::twoPi);

        Q = tensor
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
    }

    molCloud_.createMolecule
    (
        position,
        cell,   
        Q,
        v,
        vector::zero,
        pi,
        vector::zero,
        specialPosition,
        special,
        id,
        1.0,
        molCloud_.getTrackingNumber()
    );
}*/

vector polyHydroxylSiO2SurfaceBounded::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return sqrt(molCloud_.redUnits().kB()*temperature/mass)*vector
    (
        molCloud_.rndGen().GaussNormal(),
        molCloud_.rndGen().GaussNormal(),
        molCloud_.rndGen().GaussNormal()
    );
}

vector polyHydroxylSiO2SurfaceBounded::equipartitionAngularMomentum
(
    scalar temperature,
    label id
)
{
    const polyMolecule::constantProperties& cP = molCloud_.constProps(id);

    scalar sqrtKbT = sqrt(molCloud_.redUnits().kB()*temperature);

    if (cP.linearMolecule())
    {
        return sqrtKbT*vector
        (
            0.0,
            sqrt(cP.momentOfInertia().yy())*molCloud_.rndGen().GaussNormal(),
            sqrt(cP.momentOfInertia().zz())*molCloud_.rndGen().GaussNormal()
        );
    }
    else
    {
        return sqrtKbT*vector
        (
            sqrt(cP.momentOfInertia().xx())*molCloud_.rndGen().GaussNormal(),
            sqrt(cP.momentOfInertia().yy())*molCloud_.rndGen().GaussNormal(),
            sqrt(cP.momentOfInertia().zz())*molCloud_.rndGen().GaussNormal()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyHydroxylSiO2SurfaceBounded::~polyHydroxylSiO2SurfaceBounded()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



} // End namespace Foam

// ************************************************************************* //
