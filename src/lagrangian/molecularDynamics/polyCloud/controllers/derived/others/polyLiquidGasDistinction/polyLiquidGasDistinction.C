/*---------------------------------------------------------------------------*\
 *  =========                 |
 *  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 *   \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
 *     \\/     M anipulation  |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 * 
 *    OpenFOAM is free software; you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published by the
 *    Free Software Foundation; either version 2 of the License, or (at your
 *    option) any later version.
 * 
 *    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *    for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with OpenFOAM; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * Description
 * 
 * \*---------------------------------------------------------------------------*/

#include "polyLiquidGasDistinction.H"
#include "addToRunTimeSelectionTable.H"
//#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    
defineTypeNameAndDebug(polyLiquidGasDistinction, 0);
addToRunTimeSelectionTable(polyStateController, polyLiquidGasDistinction, dictionary);

void polyLiquidGasDistinction::evaluatePair
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
    label idI = molI->id();
    label idJ = molJ->id();
        
    if((idI == molIdA_ || idI == molIdB_) && (idJ == molIdA_ || idJ == molIdB_))
    {   
        vector rsIsJ = molI->position() - molJ->position();
        
        scalar rsIsJMag = mag(rsIsJ);
        
        if(rsIsJMag < rMax_)
        {
            lg_[molI->trackingNumber()]++;
            lg_[molJ->trackingNumber()]++;
        }
    }
}


void polyLiquidGasDistinction::evaluatePairS
(
    polyMolecule* molReal,
    polyMolecule* molRef
)
{
    label idI = molReal->id();
    label idJ = molRef->id();
    if((idI == molIdA_ || idI == molIdB_) && (idJ == molIdA_ || idJ == molIdB_))
    {   
        vector rsRealsRef =molReal->position() - molRef->position();            
        
        scalar rsRealsRefMag = mag(rsRealsRef);
        
        if(rsRealsRefMag < rMax_)
        {
            lg_[molReal->trackingNumber()]++;
        }
    }
}

void polyLiquidGasDistinction::readProperties()
{
    liquidLimit_ = readLabel(propsDict_.lookup("NumberOfAtomsForLiquid"));
    rMax_ = readScalar(propsDict_.lookup("rMax"));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyLiquidGasDistinction::polyLiquidGasDistinction
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    rMax_(readScalar(propsDict_.lookup("rMax"))),
    liquidLimit_(readLabel(propsDict_.lookup("NumberOfAtomsForLiquid"))),
    molIdA_(-1),
    molIdB_(-1),
    lg_()
{
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);
    
    if(regionId_ == -1)
    {
        FatalErrorIn("polyLiquidGasDistinction::polyLiquidGasDistinction()")
        << "Cannot find region: " << regionName_ << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
    
    const List<word>& idList(molCloud_.cP().molIds());
    
    const word molIdA = propsDict_.lookup("molIdA");
    molIdA_ = findIndex(idList, molIdA);

    const word molIdB = propsDict_.lookup("molIdB");
    molIdB_ = findIndex(idList, molIdB);
    
    if(molIdA_ == -1)
    {
        FatalErrorIn("polyLiquidGasDistinction::polyLiquidGasDistinction()")
        << "Cannot find molIdA: " << molIdA << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
    
    if(molIdB_ == -1)
    {
        FatalErrorIn("polyLiquidGasDistinction::polyLiquidGasDistinction()")
        << "Cannot find molIdB: " << molIdB << nl << "in: "
        << time_.time().system()/"fieldPropertiesDict"
        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyLiquidGasDistinction::~polyLiquidGasDistinction()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyLiquidGasDistinction::initialConfiguration()
{
	label lg_size = molCloud_.size();
    
    if (Pstream::parRun())
    {
        
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << lg_size;
                }
            }
        }
        
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label lgProc;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> lgProc;
                }
                
                lg_size += lgProc;
            }
        }
    }
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(mol().trackingNumber() > lg_size)
        {                
            lg_size=mol().trackingNumber();
        }
    }
    
    if (Pstream::parRun())
    {
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << lg_size;
                }
            }
        }
        
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label lgProc;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> lgProc;
                }
                
                if(lgProc > lg_size)
                    lg_size = lgProc;
            }
        }
        
    }
    
    lg_.setSize(lg_size+1,0);
    
}

void polyLiquidGasDistinction::controlBeforeVelocityI()
{}

void polyLiquidGasDistinction::controlBeforeMove()
{}


void polyLiquidGasDistinction::controlBeforeForces()
{}

void polyLiquidGasDistinction::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyLiquidGasDistinction::controlAfterForces()
{}



void polyLiquidGasDistinction::controlAfterVelocityII()
{
    Info << "polyLiquidGasDistinction: control" << endl;
    
    if(control_)
    {
        // Set the initial groupId for all of atoms in the defined zoneName and molId.
        const List< DynamicList<polyMolecule*> >& cellOccupancy = molCloud_.cellOccupancy();
        
        const labelList& cells = mesh_.cellZones()[regionId_];
        
        forAll(cells, c)
        {
            const label& cellI = cells[c];
            
            const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];
            
            forAll(molsInCell, mIC)
            {
                polyMolecule* molI = molsInCell[mIC];
                
                if(molI->id() == molIdA_)
                {
                    // Change all the liquid particles back to gas
                    molI->id() = molIdB_; 
                }
            }
        }
        
        forAll(lg_, l)
        {
            lg_[l] = 0;
        }
        
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        {
            
            polyMolecule* molI = &mol();
            
            polyMolecule* molJ = &mol();
            
            const labelListList& dil=molCloud_.il().dil();
            
            molCloud_.il().setRIPL();
            
            const List<DynamicList<polyMolecule*> >& ril=molCloud_.il().ripl();
            
            
            forAll(dil, d)
            {   
                bool flag = false;
                forAll(cells, c)
                {
                    if(cells[c] == d)
                    {
                        flag = true;
                        break;
                    }
                }
                if(flag)
                {
                    forAll(cellOccupancy[d],cellIMols)
                    {
                        molI = cellOccupancy[d][cellIMols];
                        
                        forAll(dil[d], interactingCells)
                        {
                            List< polyMolecule* > cellJ =
                            cellOccupancy[dil[d][interactingCells]];
                            
                            forAll(cellJ, cellJMols)
                            {
                                
                                molJ = cellJ[cellJMols];
                                evaluatePair(molI, molJ);
                            }
                        }
                        
                        forAll(cellOccupancy[d],cellIOtherMols)
                        {
                            molJ = cellOccupancy[d][cellIOtherMols];
                            
                            if (molJ > molI)
                            {
                                evaluatePair(molI, molJ);
                            }
                        }
                        
                        forAll(ril[d], r)
                        {
                            molJ= ril[d][r];
                            evaluatePairS(molI, molJ);
                        }
                    } // end of forAll(cell0ccupancy[d], cellImols)
                } // end of if(flag)
            }  // end of forAll(dil, d)
            
            
        } // end of iterator
        
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(mol().id() == molIdB_)
                if(lg_[mol().trackingNumber()] > liquidLimit_)
                {
                    mol().id()=molIdA_;
                }
        }
        
        Info << "PolyLiquidGasDistinction: atoms have many neighbours are identified as liquid." << endl;
    }
    
}



void polyLiquidGasDistinction::calculateProperties()
{}

void polyLiquidGasDistinction::output
(
    const fileName& fixedPathName, 
    const fileName& timePath
)
{}

void polyLiquidGasDistinction::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateStateControllerProperties(newDict);
    
    propsDict_ = newDict.subDict(typeName + "Properties");
    
    readProperties();
}
    
} // End namespace Foam

// ************************************************************************* //
