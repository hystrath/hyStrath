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
    waterMoleculeCloud
    // OBSOLETE - keeping here just in case
Description

\*----------------------------------------------------------------------------*/

#include "waterMoleculeCloud.H"
#include "potential.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- for a molecule i, perform the intermolecular force calculation, with
//  molecules within the same cell, within the real interaction cells, and referred
//  interaction cells. 
void waterMoleculeCloud::calculateForceOnMolecule
(
    waterMolecule* mol
)
{
    mol->a() = vector::zero;
    mol->potentialEnergy() = 0.0;
    mol->rf() = tensor::zero;
    mol->R() = GREAT;

    mol->siteForces() = vector::zero;

    const label& cell = mol->cell();
    const label& idI = mol->id();
    const waterMolecule::constantProperties& constPropI = constProps(idI);

//     List<label> siteIdsI = constPropI.siteIds();
//     List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();
//     List<bool> electrostaticSitesI = constPropI.electrostaticSites();    

    const pairPotentialList& pairPot(pot_.pairPotentials());
    const waterReferredCellList& referredInteractionList = il_.ril();
    const labelListList& fullInteractionCellList = il_.dil().fil();
    const pairPotential& electrostatic = pairPot.electrostatic();

    //- real cells

    label idJ = -1;

    const labelList& dICL = fullInteractionCellList[cell];

    //- molecules within direct interaction cells (not incl. owner cell)
    forAll(dICL, dCell)
    {
        const List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            const waterMolecule* molJ = molsCellJ[m];

            idJ = molJ->id();

            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();
    
                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        mol->siteForces()[sI] += fsIsJ;
    
                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
    
                        mol->rf() += virialContribution;
                    }
                }
            }

            {
                vector rIJ = mol->position() - molJ->position();
        
                scalar rIJMag = mag(rIJ);
        
                if(mol->R() > rIJMag)
                {
                    mol->R() = rIJMag;
                }
            }
        
            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        scalar chargeI = constPropI.sites()[sI].siteCharge();
        
                        scalar chargeJ = constPropJ.sites()[sJ].siteCharge();              
            
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                            *chargeI*chargeJ
                            *electrostatic.force(rsIsJMag);
            
                        mol->siteForces()[sI] += fsIsJ;
    
                        scalar potentialEnergy
                        (
                            chargeI*chargeJ
                            *electrostatic.energy(rsIsJMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
    
                        mol->rf() += virialContribution;
                    }
                }
            }
        }
    }

    //- molecules within owner cell

    forAll(cellOccupancy_[cell], mols)
    {
        const waterMolecule* molJ = cellOccupancy_[cell][mols];

        idJ = molJ->id();

        //- warning: to generalise this function, we need to test whether molJ = mol
        if(molJ != mol)
        {
            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

//             List<label> siteIdsJ = constPropJ.siteIds();
        
//             List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();
//             List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        mol->siteForces()[sI] += fsIsJ;
    
                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
    
                        mol->rf() += virialContribution;
                    }
                }
            }

            {
                vector rIJ = mol->position() - molJ->position();
        
                scalar rIJMag = mag(rIJ);
        
                if(mol->R() > rIJMag)
                {
                    mol->R() = rIJMag;
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        scalar chargeI = constPropI.sites()[sI].siteCharge();

                        scalar chargeJ = constPropJ.sites()[sJ].siteCharge();              
            
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                            *chargeI*chargeJ
                            *electrostatic.force(rsIsJMag);
            
                        mol->siteForces()[sI] += fsIsJ;
    
                        scalar potentialEnergy
                        (
                            chargeI*chargeJ
                            *electrostatic.energy(rsIsJMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
    
                        mol->rf() += virialContribution;
                    }
                }
            }
        }
    }

    //- referred cells 

    const labelList& refCellIds = referredInteractionList.refCellIds()[cell];

    forAll(refCellIds, r)
    {
        const label& refCellId = refCellIds[r];
        const waterReferredCell& refCellI = referredInteractionList[refCellId];

        forAll(refCellI, refMols)
        {
            const waterReferredMolecule* molRef = &(refCellI[refMols]);

            label idRef = molRef->id();
            const waterMolecule::constantProperties& constPropRef(constProps(idRef));
//             List<label> siteIdsRef = constPropRef.siteIds();
//             List<bool> pairPotentialSitesRef = constPropRef.pairPotentialSites();
//             List<bool> electrostaticSitesRef = constPropRef.electrostaticSites();

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropRef.pairPotSites(), pJ)
                {
                    label sRef = constPropRef.pairPotSites()[pJ];

                    vector rsRealsRef =
                        mol->sitePositions()[sI]
                    - molRef->sitePositions()[sRef];
    
                    scalar rsRealsRefMagSq = magSqr(rsRealsRef);
        
                    label idsI = constPropI.sites()[sI].siteId();
                    label idsRef = constPropRef.sites()[sRef].siteId();

                    if (pairPot.rCutSqr(idsI, idsRef, rsRealsRefMagSq))
                    {
                        scalar rsRealsRefMag = mag(rsRealsRef);
    
                        vector fsRealsRef =
                            (rsRealsRef/rsRealsRefMag)
                        *pairPot.force(idsI, idsRef, rsRealsRefMag);
    
                        mol->siteForces()[sI] += fsRealsRef;
    
                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsRef, rsRealsRefMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rRealRef = mol->position() - molRef->position();
    
                        mol->rf() +=
                            (rsRealsRef*fsRealsRef)
                        *(rsRealsRef & rRealRef)
                        /rsRealsRefMagSq;
                    }
                }
            }

            {
                vector rIJ = mol->position() - molRef->position();
        
                scalar rIJMag = mag(rIJ);
        
                if(mol->R() > rIJMag)
                {
                    mol->R() = rIJMag;
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropRef.electrostaticSites(), pJ)
                {
                    label sRef = constPropRef.electrostaticSites()[pJ];

                    vector rsRealsRef =
                        mol->sitePositions()[sI]
                        - molRef->sitePositions()[sRef];
    
                    scalar rsRealsRefMagSq = magSqr(rsRealsRef);
    
                    if (rsRealsRefMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsRealsRefMag = mag(rsRealsRef);
    
                        scalar chargeReal = constPropI.sites()[sI].siteCharge();
        
                        scalar chargeRef = constPropRef.sites()[sRef].siteCharge();                
            
                        vector fsRealsRef =
                            (rsRealsRef/rsRealsRefMag)
                            *chargeReal*chargeRef
                            *electrostatic.force(rsRealsRefMag);
    
                        mol->siteForces()[sI] += fsRealsRef;
    
                        scalar potentialEnergy
                        (
                            chargeReal*chargeRef
                            *electrostatic.energy(rsRealsRefMag)
                        );
    
                        mol->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rRealRef = mol->position() - molRef->position();
    
                        mol->rf() +=
                            (rsRealsRef*fsRealsRef)
                        *(rsRealsRef & rRealRef)
                        /rsRealsRefMagSq;
                    }
                }
            }
        }
    }
    
    mol->setSitePositions(constPropI);

    const List<vector>& siteForces = mol->siteForces();

    scalar massI = constPropI.mass();

//     mol->a() = vector::zero;
    mol->tau() = vector::zero;

    forAll(siteForces, s)
    {
        const vector& f = siteForces[s];

        mol->a() += f/massI;

        mol->tau() += (constPropI.sites()[s].siteReferencePosition() ^ (mol->Q().T() & f));
    }

//     Info << "pE: " << mol->potentialEnergy() << endl;
}


void waterMoleculeCloud::calculateShortestRadius
(
    waterMolecule* mol
//     const potential& pot
)
{

    mol->R() = GREAT;

    const label& cell = mol->cell();
    const label& idI = mol->id();
    const waterMolecule::constantProperties& constPropI = constProps(idI);

/*    List<label> siteIdsI = constPropI.siteIds();
    List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();
    List<bool> electrostaticSitesI = constPropI.electrostaticSites();  */ 

    const pairPotentialList& pairPot(pot_.pairPotentials());
    const pairPotential& electrostatic = pairPot.electrostatic();
    const waterReferredCellList& referredInteractionList = il_.ril();
    const waterDirectInteractionList& directInteractionCellList = il_.dil();
    const labelListList& fullInteractionCellList = directInteractionCellList.fil();


    //- real cells

    label idJ = -1;

    const labelList& dICL = fullInteractionCellList[cell];

    //- molecules within direct interaction cells (not incl. owner cell)
    forAll(dICL, dCell)
    {
        const List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            const waterMolecule* molJ = molsCellJ[m];

            idJ = molJ->id();

            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];
        
                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        if(mol->R() > rsIsJMag)
                        {
                            mol->R() = rsIsJMag;
                        }
                    }
                }
            }


            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        if(mol->R() > rsIsJMag)
                        {
                            mol->R() = rsIsJMag;
                        }
                    }
                }
            }
        }
    }
      
    //- molecules within owner cell
    forAll(cellOccupancy_[cell], mols)
    {
        const waterMolecule* molJ = cellOccupancy_[cell][mols];

        idJ = molJ->id();

        //- warning: to generalise this function, we need to test whether molJ = mol
        if(molJ != mol)
        {
            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        if(mol->R() > rsIsJMag)
                        {
                            mol->R() = rsIsJMag;
                        }
                    }
                }
            }

//                     if (electrostaticSitesI[sI] && electrostaticSitesJ[sJ])
//                     {
            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        if(mol->R() > rsIsJMag)
                        {
                            mol->R() = rsIsJMag;
                        }
                    }
                }
            }
        }
    }


    //- referred cells 

    const labelList& refCellIds = referredInteractionList.refCellIds()[cell];

    forAll(refCellIds, r)
    {
        const label& refCellId = refCellIds[r];
        const waterReferredCell& refCellI = referredInteractionList[refCellId];

        forAll(refCellI, refMols)
        {
            const waterReferredMolecule* molRef = &(refCellI[refMols]);

            label idRef = molRef->id();
            const waterMolecule::constantProperties& constPropRef(constProps(idRef));
//             List<label> siteIdsRef = constPropRef.siteIds();
//             List<bool> pairPotentialSitesRef = constPropRef.pairPotentialSites();
//             List<bool> electrostaticSitesRef = constPropRef.electrostaticSites();

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropRef.pairPotSites(), pJ)
                {
                    label sRef = constPropRef.pairPotSites()[pJ];
        
                    vector rsRealsRef =
                        mol->sitePositions()[sI]
                    - molRef->sitePositions()[sRef];
    
                    scalar rsRealsRefMagSq = magSqr(rsRealsRef);

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsRef = constPropRef.sites()[sRef].siteId();

                    if (pairPot.rCutSqr(idsI, idsRef, rsRealsRefMagSq))
                    {
                        scalar rsRealsRefMag = mag(rsRealsRef);
    
                        if(mol->R() > rsRealsRefMag)
                        {
                            mol->R() = rsRealsRefMag;
                        }
                    }
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropRef.electrostaticSites(), pJ)
                {
                    label sRef = constPropRef.electrostaticSites()[pJ];

                    vector rsRealsRef =
                        mol->sitePositions()[sI]
                    - molRef->sitePositions()[sRef];
    
                    scalar rsRealsRefMagSq = magSqr(rsRealsRef);
    
                    if (rsRealsRefMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsRealsRefMag = mag(rsRealsRef);
    
                        if(mol->R() > rsRealsRefMag)
                        {
                            mol->R() = rsRealsRefMag;
                        }
                    }
                }
            }
        }
    }
}

//- MOLECULAR INSERTION OPERATIONS


// update the forces and energies on those molecules occupying real cells,
// only within interaction range -- i.e. direct interaction list -- 
// due to the insertion of a molecule. 
void waterMoleculeCloud::updateForceOnRealCellMoleculesDueToInsertion
(
    waterMolecule* mol
//     const potential& pot
)
{

    const label& cell = mol->cell();
    const label& idI = mol->id();
    const waterMolecule::constantProperties& constPropI = constProps(idI);

//     List<label> siteIdsI = constPropI.siteIds();
//     List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();
//     List<bool> electrostaticSitesI = constPropI.electrostaticSites();

    const pairPotentialList& pairPot(pot_.pairPotentials());
//     const waterReferredCellList& referredInteractionList = il_.ril();
    const pairPotential& electrostatic = pairPot.electrostatic();

    const waterDirectInteractionList& directInteractionCellList = il_.dil();
    const labelListList& fullInteractionCellList = directInteractionCellList.fil();


    //- real cells

    label idJ = -1;

    const labelList& dICL = fullInteractionCellList[cell];

    //- molecules within direct interaction cells (not incl. owner cell)
    forAll(dICL, dCell)
    {
        const List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            waterMolecule* molJ = molsCellJ[m];

            idJ = molJ->id();

            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

            scalar massJ = constPropJ.mass();
        
            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];
        
                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
        
                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        molJ->siteForces()[sJ] += -fsIsJ;
    
                        molJ->a() += -fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsIsJ));

                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        molJ->potentialEnergy() += 0.5*potentialEnergy;

                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

    //                             tensor virialContribution = (rsIsJ*fsIsJ);
    
                        molJ->rf() += virialContribution;
                    }
                }
            }

            {
                vector rIJ = mol->position() - molJ->position();
        
                scalar rIJMag = mag(rIJ);
        
                if(molJ->R() > rIJMag)
                {
                    molJ->R() = rIJMag;
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
        
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
                    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                            *chargeI*chargeJ*electrostatic.force(rsIsJMag);
    
                        molJ->siteForces()[sJ] += -fsIsJ;
    
                        molJ->a() += -fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsIsJ));

                        scalar potentialEnergy
                        (
                                chargeI*chargeJ
                                *electrostatic.energy(rsIsJMag)
                        );
    
                        molJ->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

//                             tensor virialContribution = (rsIsJ*fsIsJ);

                        molJ->rf() += virialContribution;
                    }
                }
            }
        }
    }
      
    //- molecules within owner cell
    forAll(cellOccupancy_[cell], mols)
    {
        waterMolecule* molJ = cellOccupancy_[cell][mols];

        idJ = molJ->id();

        //- warning: to generalise this function, we need to test whether molJ = mol
        if(molJ != mol)
        {
            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

//             List<label> siteIdsJ = constPropJ.siteIds();
//         
//             List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();
//             List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

            scalar massJ = constPropJ.mass();
        
//             forAll(siteIdsI, sI)
//             {
//                 label idsI(siteIdsI[sI]);
//         
//                 forAll(siteIdsJ, sJ)
//                 {
//                     label idsJ(siteIdsJ[sJ]);

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
        
                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        molJ->siteForces()[sJ] += -fsIsJ;
    
                        molJ->a() += -fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsIsJ));

                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        molJ->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

    //                             tensor virialContribution = (rsIsJ*fsIsJ);
    
                        molJ->rf() += virialContribution;
                    }
                }
            }

            {
                vector rIJ = mol->position() - molJ->position();
        
                scalar rIJMag = mag(rIJ);
        
                if(molJ->R() > rIJMag)
                {
                    molJ->R() = rIJMag;
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    scalar chargeI = constPropI.sites()[sI].siteCharge();

                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
                    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                            *chargeI*chargeJ*electrostatic.force(rsIsJMag);
    
                        molJ->siteForces()[sJ] += -fsIsJ;
    
                        molJ->a() += -fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsIsJ));

                        scalar potentialEnergy
                        (
                                chargeI*chargeJ
                                *electrostatic.energy(rsIsJMag)
                        );
    
                        molJ->potentialEnergy() += 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

//                             tensor virialContribution = (rsIsJ*fsIsJ);

                        molJ->rf() += virialContribution;
                    }
                }
            }
        }
    }
}


// update the forces and energies on those molecules occupying referred cells,
// only within interaction range due to the insertion of a molecule. 
// Here we only consider referred cells belonging to periodic patches, since
// the molecules are available on the same mesh and hence no parallel processing
// is necessary. Does not consider those referred cells belonging to processor patches
// void moleculeCloud::updateForceOnReferredCellMoleculesDueToInsertion
// (
//     molecule* mol
// )
// {
//     const label& cell = mol->cell();
//     const label& idI = mol->id();
//     const vector& posI = mol->position();
// 
//     vector rIJ;
//     scalar rIJMag;
//     scalar rIJMagSq;
//     scalar fIJMag;
//     label idJ;
// //     scalar massJ;
// 
//     //- referred cells 
// 
//     const labelList& refCellIds = referredInteractionList_.refCellIds()[cell];
// 
//     forAll(refCellIds, r)
//     {
//         const label& refCellId = refCellIds[r];
//         waterReferredCell& refCellI = referredInteractionList_[refCellId];
// 
//         if(refCellI.sourceProc() == Pstream::myProcNo())
//         {
//             const label& sourceCell = refCellI.sourceCell();
// 
//             List<molecule*>& molsInCell = cellOccupancy_[sourceCell];
// 
//             forAll(refCellI, r)
//             {
//                 waterReferredMolecule* molJ = &(refCellI[r]);
//             
//                 molecule* molR = molsInCell[r];
//     
//                 idJ = molJ->id();
//                 rIJ = posI - molJ->position();
//                 
//                 rIJMagSq = magSqr(rIJ);
//                 
//                 if(pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
//                 {
//                     rIJMag = mag(rIJ);
//                     
//                     fIJMag = pairPotentials_.force(idI, idJ, rIJMag);
//                     
//                     // Acceleration increment for mol
//                     molR->A() += rIJ*fIJMag/(molR->mass());
//                     
//                     scalar potentialEnergy
//                     (
//                         pairPotentials_.energy(idI, idJ, rIJMag)
//                     );
//                     
//                     molR->pE() += 0.5*potentialEnergy;
//                     
//                     molR->rDotf() += 0.5*fIJMag*rIJMagSq;
//                     
//                     //- shortest radius between molecules
//                     if(molR->R() > rIJMag)
//                     {
//                         molR->R() = rIJMag;
//                     }
//                 }
//             }
//         }
//     }
// }


//- insert molecule within the referred cells and update the forces of
//  molecules within the list of real interaction cells 
//  if procNo is myProcNo(), then it will treat referred cells associated
//  with periodic patches.
void waterMoleculeCloud::updateForceOnReferredCellMoleculesDueToInsertion
(
    const label& cellId,
    const label& idI,
    const vector& posI,
    const List<vector>& sitePos,
    const label& cellRecId
//     const potential& pot
)
{
//     const waterMolecule::constantProperties& constPropI = constProps(idI);


    const pairPotentialList& pairPot(pot_.pairPotentials());
    const pairPotential& electrostatic = pairPot.electrostatic(); 

    waterReferredCellList& referredInteractionList = il_.ril();
    const List<receivingReferralList>& cellReceivingReferralLists = il_.cellReceivingReferralLists();

    const labelList& refCellsId = cellReceivingReferralLists[cellRecId][cellId];

    label idJ = -1;

    forAll(refCellsId, r)
    {
        const label& refCellId = refCellsId[r];

        waterReferredCell& refCellI = referredInteractionList[refCellId];//*****

        refCellI.referInMolecule
        (
            waterReferredMolecule(idI, posI, sitePos, vector::zero, 1.0, -1)
        );

        //- update forces on real molecules within the cells of interaction

        waterReferredMolecule* molR = &(refCellI[refCellI.size()-1]);

        label idR = molR->id();

        const waterMolecule::constantProperties& constPropR = constProps(idR);
   

        const labelList& realCells = refCellI.realCellsForInteraction();
        
        forAll(realCells, rC)
        {
            const label& realCellI = realCells[rC];

            List<waterMolecule*>& molsInCell = cellOccupancy_[realCellI];

            forAll(molsInCell, m)
            {
                waterMolecule* molJ = molsInCell[m];
    
                idJ = molJ->id();

                const waterMolecule::constantProperties& constPropJ(constProps(idJ));
    
                scalar massJ = constPropJ.mass();
            
                forAll(constPropR.pairPotSites(), pI)
                {
                    label sR = constPropR.pairPotSites()[pI];
            
                    forAll(constPropJ.pairPotSites(), pJ)
                    {
                        label sJ = constPropJ.pairPotSites()[pJ];

                        label idsR = constPropR.sites()[sR].siteId();
                        label idsJ = constPropJ.sites()[sJ].siteId();

                        vector rsRsJ =
                            molR->sitePositions()[sR] - molJ->sitePositions()[sJ];
        
                        scalar rsRsJMagSq = magSqr(rsRsJ);
        
                        if(pairPot.rCutSqr(idsR, idsJ, rsRsJMagSq))
                        {
                            scalar rsRsJMag = mag(rsRsJ);
        
                            vector fsRsJ =
                                (rsRsJ/rsRsJMag)
                            *pairPot.force(idsR, idsJ, rsRsJMag);

                            // - update site forces            
                            molJ->siteForces()[sJ] += -fsRsJ;
        
                            // - update whole waterMolecule

                            molJ->a() += -fsRsJ/massJ;
    
                            molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsRsJ));

                            scalar potentialEnergy
                            (
                                pairPot.energy(idsR, idsJ, rsRsJMag)
                            );

                            molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                            vector rRJ = molR->position() - molJ->position();
        
                            tensor virialContribution =
                                (rsRsJ*fsRsJ)*(rsRsJ & rRJ)/rsRsJMagSq;

//                                 tensor virialContribution = (rsRsJ*fsRsJ);

                            molJ->rf() += virialContribution;
                        }
                    }
                }

                {
                    vector rIJ = molR->position() - molJ->position();
            
                    scalar rIJMag = mag(rIJ);
            
                    if(molJ->R() > rIJMag)
                    {
                        molJ->R() = rIJMag;
                    }
                }

                forAll(constPropR.electrostaticSites(), pI)
                {
                    label sR = constPropR.electrostaticSites()[pI];
            
                    forAll(constPropJ.electrostaticSites(), pJ)
                    {
                        label sJ = constPropJ.electrostaticSites()[pJ];

                        vector rsRsJ =
                            molR->sitePositions()[sR] - molJ->sitePositions()[sJ];
        
                        scalar rsRsJMagSq = magSqr(rsRsJ);
        
                        if(rsRsJMagSq <= electrostatic.rCutSqr())
                        {
                            scalar rsRsJMag = mag(rsRsJ);

                            scalar chargeR = constPropR.sites()[sR].siteCharge();
                            scalar chargeJ = constPropJ.sites()[sJ].siteCharge();

                            vector fsRsJ =
                                (rsRsJ/rsRsJMag)
                                *chargeR*chargeJ*electrostatic.force(rsRsJMag);

                            // - update site forces            
                            molJ->siteForces()[sJ] += -fsRsJ;
        
                            // - update whole molecule

                            molJ->a() += -fsRsJ/massJ;
    
                            molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & -fsRsJ));

                            scalar potentialEnergy
                            (
                                chargeR*chargeJ
                                *electrostatic.force(rsRsJMag)
                            );

                            molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                            vector rRJ = molR->position() - molJ->position();
        
                            tensor virialContribution =
                                (rsRsJ*fsRsJ)*(rsRsJ & rRJ)/rsRsJMagSq;
//                                 tensor virialContribution = rsRsJ*fsRsJ;

                            molJ->rf() += virialContribution;
                        }
                    }
                }
            }
        }
    }
}





// update the forces and energies on those molecules occupying referred cells,
// only within interaction range due to the insertion of a molecule. 
// Here we are only considering referred cells from processor patches.
// This function must be called by all processors due to communication purposes.
// void moleculeCloud::updateForceOnRealCellMoleculesDueToInsertionParallel
// (
//     const label& cellI,
//     const label& procN
// )
// {
//     List< DynamicList<label> > cellsOnProcs(Pstream::nProcs());
// 
//     if(procN == Pstream::myProcNo())
//     {
//         molecule* newMol = last();
//         const label& molId = newMol->id();
//         const vector& molPosition = newMol->position();
// 
//         const labelList& cellRefIds = referredInteractionList_.refCellIds()[cellI];
//        
//         forAll(cellRefIds, c)
//         {
//             const waterReferredCell& refCellI = referredInteractionList_[cellRefIds[c]];
// 
//             if(refCellI.sourceProc() != Pstream::myProcNo())
//             {
//                 cellsOnProcs[refCellI.sourceProc()].append(refCellI.sourceCell());
//             }
//         }
// 
//         //- sending
//         for (int i = 0; i < Pstream::nProcs(); i++)
//         {
//             if(i != Pstream::myProcNo())
//             {
//                 const int proc = i;
//                 {
//                     OPstream toNeighbour(proc);
//                     toNeighbour << cellsOnProcs[i].shrink() << molId << molPosition << cellI;
//                 }
//             }
//         }
//     }
// 
//     //- receiving
//     if(procN != Pstream::myProcNo())
//     {
//         List<label> cellLabels;
//         label idI;
//         vector posI;
//         label cellFromProc;
// 
//         const int proc = procN;
//         {
//             IPstream fromNeighbour(proc);
//             fromNeighbour >> cellLabels >> idI >> posI >> cellFromProc;
//         }
// 
//         vector rIJ;
//         scalar rIJMag;
//         scalar rIJMagSq;
//         scalar fIJMag;
//         label idJ;
// //         scalar massJ;
// 
//         forAll(cellLabels, c)
//         {
//             const label& cell = cellLabels[c];
// 
//             forAll(cellOccupancy_[cell], mols)
//             {
//                 molecule* molJ = cellOccupancy_[cell][mols];
//         
//                 idJ = molJ->id();
// //                 massJ = molJ->mass();
//                 rIJ = posI - molJ->position();
//                 
//                 rIJMagSq = magSqr(rIJ);
//                 
//                 if(pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
//                 {
//                     rIJMag = mag(rIJ);
//                 
//                     fIJMag = pairPotentials_.force(idI, idJ, rIJMag);
//                 
//                     // Acceleration increment for mol
//                     molJ->A() += rIJ*fIJMag/(molJ->mass());
//                 
//                     scalar potentialEnergy
//                     (
//                         pairPotentials_.energy(idI, idJ, rIJMag)
//                     );
//                 
//                     molJ->pE() += 0.5*potentialEnergy;
//                 
//                     molJ->rDotf() += 0.5*fIJMag*rIJMagSq;
//                 
//                     //- shortest radius between molecules
//                     if(molJ->R() > rIJMag)
//                     {
//                         molJ->R() = rIJMag;
//                     }
//                 }
//             }
//         }
// 
//         //- place molecule in referred cell
//         const labelList& cellRefIds = referredInteractionList_.refCellIds()[cellLabels[0]];
//         
//         forAll(cellRefIds, c)
//         {
//             waterReferredCell& refCellI = referredInteractionList_[cellRefIds[c]];
// 
//             if((refCellI.sourceProc() == procN) && (refCellI.sourceCell() == cellFromProc))
//             {
//                 refCellI.append
//                 (
//                     waterReferredMolecule(idI, posI)
//                 );
// 
//                 refCellI.shrink();
//             }
//         }
//     }
// }


//- MOLECULAR DELETION OPERATIONS

// update the forces and energies on those molecules occupying real cells,
// only within interaction range -- i.e. direct interaction list -- 
// due to the deletion of a molecule. 
void waterMoleculeCloud::updateForceOnRealCellMoleculesDueToDeletion
(
    waterMolecule* mol
)
{
	//- here we assume that the empty occupancy volume of the
	// neighbouring molecules is re-built in the next time-step.

    const label& cell = mol->cell();
    const label& idI = mol->id();
    const waterMolecule::constantProperties& constPropI = constProps(idI);

//     List<label> siteIdsI = constPropI.siteIds();
//     List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();
//     List<bool> electrostaticSitesI = constPropI.electrostaticSites();

    const pairPotentialList& pairPot(pot_.pairPotentials());
    const pairPotential& electrostatic = pairPot.electrostatic();    

//     const waterReferredCellList& referredInteractionList = il_.ril();
    const waterDirectInteractionList& directInteractionCellList = il_.dil();
    const labelListList& fullInteractionCellList = directInteractionCellList.fil();


    //- real cells

    label idJ = -1;

    const labelList& dICL = fullInteractionCellList[cell];

    //- molecules within direct interaction cells (not incl. owner cell)
    forAll(dICL, dCell)
    {
        const List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            waterMolecule* molJ = molsCellJ[m];

            idJ = molJ->id();

            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

//             List<label> siteIdsJ = constPropJ.siteIds();
//         
//             List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();
//             List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

            scalar massJ = constPropJ.mass();
        
            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        molJ->siteForces()[sJ] += fsIsJ;
    
                        molJ->a() += fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsIsJ));

                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        molJ->potentialEnergy() -= 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;


//                             tensor virialContribution = (rsIsJ*fsIsJ);

                        molJ->rf() -= virialContribution;
                    }
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
        
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
                    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                        scalar rsIsJMag = mag(rsIsJ);
                
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                            *chargeI*chargeJ*electrostatic.force(rsIsJMag);
    
//                         molJ->siteForces()[sJ] += -fsIsJ;
//     
//                         molJ->a() += -fsIsJ/massJ;

                        molJ->siteForces()[sJ] += fsIsJ;
    
                        molJ->a() += fsIsJ/massJ;
            
//                         molJ->tau() += (constPropJ.siteReferencePositions()[sJ] ^ (molJ->Q().T() & -fsIsJ));

                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsIsJ));

                        scalar potentialEnergy
                        (
                                chargeI*chargeJ
                            *electrostatic.energy(rsIsJMag)
                        );
    
//                         molJ->potentialEnergy() += 0.5*potentialEnergy;
                        molJ->potentialEnergy() -= 0.5*potentialEnergy;

                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

//                             tensor virialContribution = (rsIsJ*fsIsJ);
    
                        molJ->rf() -= virialContribution;
                    }
                }
            }
        }
    }
      
    //- molecules within owner cell
    forAll(cellOccupancy_[cell], mols)
    {
        waterMolecule* molJ = cellOccupancy_[cell][mols];

        idJ = molJ->id();

        //- warning: to generalise this function, we need to test whether molJ = mol
        if(molJ != mol)
        {
            const waterMolecule::constantProperties& constPropJ(constProps(idJ));

//             List<label> siteIdsJ = constPropJ.siteIds();
//         
//             List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();
//             List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

            scalar massJ = constPropJ.mass();
        
//             forAll(siteIdsI, sI)
//             {
//                 label idsI(siteIdsI[sI]);
//         
//                 forAll(siteIdsJ, sJ)
//                 {
//                     label idsJ(siteIdsJ[sJ]);
//         
//                     if (pairPotentialSitesI[sI] && pairPotentialSitesJ[sJ])
//                     {

            forAll(constPropI.pairPotSites(), pI)
            {
                label sI = constPropI.pairPotSites()[pI];
        
                forAll(constPropJ.pairPotSites(), pJ)
                {
                    label sJ = constPropJ.pairPotSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);

                    label idsI = constPropI.sites()[sI].siteId();
                    label idsJ = constPropJ.sites()[sJ].siteId();

                    if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                    {
                        scalar rsIsJMag = mag(rsIsJ);
    
                        vector fsIsJ =
                            (rsIsJ/rsIsJMag)
                        *pairPot.force(idsI, idsJ, rsIsJMag);
    
                        molJ->siteForces()[sJ] += fsIsJ;
    
                        molJ->a() += fsIsJ/massJ;
            
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsIsJ));

                        scalar potentialEnergy
                        (
                            pairPot.energy(idsI, idsJ, rsIsJMag)
                        );
    
                        molJ->potentialEnergy() -= 0.5*potentialEnergy;
    
                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
//                             tensor virialContribution = (rsIsJ*fsIsJ);
    
                        molJ->rf() -= virialContribution;
                    }
                }
            }

            forAll(constPropI.electrostaticSites(), pI)
            {
                label sI = constPropI.electrostaticSites()[pI];
        
                forAll(constPropJ.electrostaticSites(), pJ)
                {
                    label sJ = constPropJ.electrostaticSites()[pJ];

                    vector rsIsJ =
                        mol->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                    scalar rsIsJMagSq = magSqr(rsIsJ);
    
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
        
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
                    
                    if(rsIsJMagSq <= electrostatic.rCutSqr())
                    {
                                    scalar rsIsJMag = mag(rsIsJ);
                
                        vector fsIsJ = (rsIsJ/rsIsJMag)*chargeI*chargeJ*electrostatic.force(rsIsJMag);
        
//                         molJ->siteForces()[sJ] += -fsIsJ;
//     
//                         molJ->a() += -fsIsJ/massJ;
                        molJ->siteForces()[sJ] += fsIsJ;
    
                        molJ->a() += fsIsJ/massJ;

//                         molJ->tau() += (constPropJ.siteReferencePositions()[sJ] ^ (molJ->Q().T() & -fsIsJ));
                        molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsIsJ));

                        scalar potentialEnergy
                        (
                                chargeI*chargeJ
                            *electrostatic.energy(rsIsJMag)
                        );
    
//                         molJ->potentialEnergy() += 0.5*potentialEnergy;

                        molJ->potentialEnergy() -= 0.5*potentialEnergy;

                        vector rIJ = mol->position() - molJ->position();
    
                        tensor virialContribution =
                            (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;

//                             tensor virialContribution = (rsIsJ*fsIsJ);

//                         molJ->rf() += virialContribution;
                        molJ->rf() -= virialContribution;
                    }
                }
            }
        }
    }
}



// update the forces and energies on those molecules occupying referred cells,
// only within interaction range due to the deletion of a molecule. 
// Here we are only considering referred cells from periodic patches, since these
// cells and their molecules are available on the same mesh unlike 
// processor-referred-cells
// void moleculeCloud::updateForceOnReferredCellMoleculesDueToDeletion
// (
//     molecule* mol
// )
// {
//     const label& cell = mol->cell();
//     const label& idI = mol->id();
//     const vector& posI = mol->position();
// 
//     vector rIJ;
//     scalar rIJMag;
//     scalar rIJMagSq;
//     scalar fIJMag;
//     label idJ;
// //     scalar massJ;
// 
//     //- referred cells 
// 
//     const labelList& refCellIds = referredInteractionList_.refCellIds()[cell];
// 
//     forAll(refCellIds, r)
//     {
//         const label& refCellId = refCellIds[r];
//         waterReferredCell& refCellI = referredInteractionList_[refCellId];
// 
//         if(refCellI.sourceProc() == Pstream::myProcNo())
//         {
//             const label& sourceCell = refCellI.sourceCell();
// 
//             List<molecule*>& molsInCell = cellOccupancy_[sourceCell];
// 
//             forAll(refCellI, r)
//             {
//                 waterReferredMolecule* molJ = &(refCellI[r]);
//             
//                 molecule* molR = molsInCell[r];
//     
//                 idJ = molJ->id();
// //                 massJ = molR->mass();
//                 rIJ = posI - molJ->position();
//                 
//                 rIJMagSq = magSqr(rIJ);
//                 
//                 if(pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
//                 {
//                     rIJMag = mag(rIJ);
//                     
//                     fIJMag = pairPotentials_.force(idI, idJ, rIJMag);
//                     
//                     
//                     // Acceleration increment for mol
//                     molR->A() -= rIJ*fIJMag/(molR->mass());
//                     
//                     scalar potentialEnergy
//                     (
//                         pairPotentials_.energy(idI, idJ, rIJMag)
//                     );
//                     
//                     molR->pE() -= 0.5*potentialEnergy;
//                     
//                     molR->rDotf() -= 0.5*fIJMag*rIJMagSq;
//                 }
//             }
//         }
//     }
// }

void waterMoleculeCloud::updateForceOnReferredCellMoleculesDueToDeletion
(
    const label& cellMolRemoveId,
    const label& cellId,
    const label& idI,
    const vector& posI,
    const label& cellRecId
//     const potential& pot
)
{
//     const waterMolecule::constantProperties& constPropI = constProps(idI);

//     List<label> siteIdsI = constPropI.siteIds();
//     List<bool> pairPotentialSitesI = constPropI.pairPotentialSites();
//     List<bool> electrostaticSitesI = constPropI.electrostaticSites();

    const pairPotentialList& pairPot(pot_.pairPotentials());
    const pairPotential& electrostatic = pairPot.electrostatic();

    waterReferredCellList& referredInteractionList = il_.ril();
    const List<receivingReferralList>& cellReceivingReferralLists = il_.cellReceivingReferralLists();

    const labelList& refCellsId = cellReceivingReferralLists[cellRecId][cellId];

    label idJ = -1;

    forAll(refCellsId, r)
    {
        const label& refCellId = refCellsId[r];

        waterReferredCell& refCellI = referredInteractionList[refCellId];

        label refMolID = -1;
        scalar minR = GREAT;

        forAll(refCellI, rM)
        {
            waterReferredMolecule* molR = &(refCellI[rM]);
    
            label idR = molR->id();
            const vector& posR = molR->position();

            const waterReferredCell& refCellITemp = referredInteractionList[refCellId];
            const vector refPos = refCellITemp.referPosition(posI);

            scalar deltaR = magSqr(refPos - posR);
    
            if((minR > deltaR) && (idI == idR))
            {
                refMolID = rM;
                minR = deltaR;
//                 Info << "incoming molecule: position: " << refPos << ", id: " << idI 
//                 << ", chosen molecule: position: " << posR << ", id: " << idR << endl;
            }
        }

        //-check (test mistake)
        if(refMolID != cellMolRemoveId)
        {
            waterReferredMolecule* molR = &(refCellI[cellMolRemoveId]);
            label idR = molR->id();
            const vector& posR = molR->position();

            const waterReferredCell& refCellITemp = referredInteractionList[refCellId];
            const vector refPos = refCellITemp.referPosition(posI);
            const label& sCI = refCellI.sourceCell();

            waterReferredMolecule* molRTrue = &(refCellI[refMolID]);
            label idRTrue = molRTrue->id();
            const vector& posRTrue = molRTrue->position();


            Info << "WARNING. Incoming molecule: position: " << refPos << ", id: " << idI 
                 << " properly chosen molecule: position " << posRTrue << ", id: " 
                 << idRTrue << ", refMolID: " << refMolID
                 << ", old chosen molecule: position: " << posR << ", id: " << idR  
                 << ", refMolID: " << cellMolRemoveId << ", waterReferredCell: " << sCI 
                 << endl;
        }


        //- update forces on real molecules within the cells of interaction
        if(refMolID != -1)
        {
            waterReferredMolecule* molR = &(refCellI[refMolID]);
    
            label idR = molR->id();
    //         const vector& posR = molR->position();
    
            const waterMolecule::constantProperties& constPropR = constProps(idR);
        
//             List<label> siteIdsR = constPropR.siteIds();
    
            const labelList& realCells = refCellI.realCellsForInteraction();
            
            forAll(realCells, rC)
            {
                const label& realCellI = realCells[rC];
    
                List<waterMolecule*>& molsInCell = cellOccupancy_[realCellI];
    
                forAll(molsInCell, m)
                {
                    waterMolecule* molJ = molsInCell[m];
        
                    idJ = molJ->id();
    //                 rIJ = posR - molJ->position();
    
                    const waterMolecule::constantProperties& constPropJ(constProps(idJ));
        
//                     List<label> siteIdsJ = constPropJ.siteIds();
//                 
//                     List<bool> pairPotentialSitesJ = constPropJ.pairPotentialSites();
//                     List<bool> electrostaticSitesJ = constPropJ.electrostaticSites();

                    scalar massJ = constPropJ.mass();

                    forAll(constPropR.pairPotSites(), pI)
                    {
                        label sR = constPropR.pairPotSites()[pI];
                
                        forAll(constPropJ.pairPotSites(), pJ)
                        {
                            label sJ = constPropJ.pairPotSites()[pJ];
                
                            vector rsRsJ =
                                molR->sitePositions()[sR] - molJ->sitePositions()[sJ];

                            label idsR = constPropR.sites()[sR].siteId();
                            label idsJ = constPropJ.sites()[sJ].siteId();         

                            scalar rsRsJMagSq = magSqr(rsRsJ);
            
                            if(pairPot.rCutSqr(idsR, idsJ, rsRsJMagSq))
                            {
                                scalar rsRsJMag = mag(rsRsJ);
            
                                vector fsRsJ =
                                    (rsRsJ/rsRsJMag)
                                *pairPot.force(idsR, idsJ, rsRsJMag);

                                // - update site forces            
                                molJ->siteForces()[sJ] += fsRsJ;
            
                                // - update whole molecule

                                molJ->a() += fsRsJ/massJ;
        
                                molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsRsJ));

                                scalar potentialEnergy
                                (
                                    pairPot.energy(idsR, idsJ, rsRsJMag)
                                );

                                molJ->potentialEnergy() -= 0.5*potentialEnergy;
            
                                vector rRJ = molR->position() - molJ->position();
            
                                tensor virialContribution =
                                    (rsRsJ*fsRsJ)*(rsRsJ & rRJ)/rsRsJMagSq;

//                                     tensor virialContribution = (rsRsJ*fsRsJ);

                                molJ->rf() -= virialContribution;

                            }
                        }
                    }
                            
                    forAll(constPropR.electrostaticSites(), pI)
                    {
                        label sR = constPropR.electrostaticSites()[pI];
                
                        forAll(constPropJ.electrostaticSites(), pJ)
                        {
                            label sJ = constPropJ.electrostaticSites()[pJ];

                            vector rsRsJ =
                                molR->sitePositions()[sR] - molJ->sitePositions()[sJ];
                    
                            scalar rsRsJMagSq = magSqr(rsRsJ);
                    
                            if(rsRsJMagSq <= electrostatic.rCutSqr())
                            {
                                scalar rsRsJMag = mag(rsRsJ);
                    
                                scalar chargeR = constPropR.sites()[sR].siteCharge();
            
                                scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
                                
                                vector fsRsJ =
                                    (rsRsJ/rsRsJMag)
                                    *chargeR*chargeJ*electrostatic.force(rsRsJMag);
            
                                // - update site forces            
//                                 molJ->siteForces()[sJ] += -fsRsJ;
                                molJ->siteForces()[sJ] += fsRsJ;

                                // - update whole molecule
            
//                                 molJ->a() += -fsRsJ/massJ;
                                molJ->a() += fsRsJ/massJ;

//                                 molJ->tau() += (constPropJ.siteReferencePositions()[sJ] ^ (molJ->Q().T() & -fsRsJ));
                                molJ->tau() += (constPropJ.sites()[sJ].siteReferencePosition() ^ (molJ->Q().T() & fsRsJ));

                                scalar potentialEnergy
                                (
                                    chargeR*chargeJ
                                    *electrostatic.force(rsRsJMag)
                                );
            
//                                 molJ->potentialEnergy() += 0.5*potentialEnergy;

                                molJ->potentialEnergy() -= 0.5*potentialEnergy;

                                vector rRJ = molR->position() - molJ->position();
                    
                                tensor virialContribution =
                                    (rsRsJ*fsRsJ)*(rsRsJ & rRJ)/rsRsJMagSq;

//                                     tensor virialContribution = (rsRsJ*fsRsJ);

//                                 molJ->rf() += virialContribution;
                                molJ->rf() -= virialContribution;
                            }
                        }
                    }
                }
            }
    
            //- remove molecule from referred cell
            refCellI.deleteReferredMolecule(refMolID);
        }
    }
}


// update the forces and energies on those molecules occupying referred cells,
// only within interaction range due to the deletion of a molecule. 
// Here we are only considering referred cells from processor patches.
// This function must be called by all processors due to communication purposes.
// void moleculeCloud::updateForceOnReferredCellMoleculesDueToDeletionParallel
// (
//     const label& cellId,
//     const label& cellI,
//     const label& procN
// )
// {
//     List< DynamicList<label> > cellsOnProcs(Pstream::nProcs());
// 
//     if(procN == Pstream::myProcNo())
//     {
//         molecule* delMol = cellOccupancy_[cellI][cellId];
//         const label& molId = delMol->id();
//         const vector& molPosition = delMol->position();
// 
//         const labelList& cellRefIds = referredInteractionList_.refCellIds()[cellI];
//        
//         forAll(cellRefIds, c)
//         {
//             const waterReferredCell& refCellI = referredInteractionList_[cellRefIds[c]];
// 
//             if(refCellI.sourceProc() != Pstream::myProcNo())
//             {
//                 cellsOnProcs[refCellI.sourceProc()].append(refCellI.sourceCell());
//             }
//         }
// 
//         //- sending
//         for (int i = 0; i < Pstream::nProcs(); i++)
//         {
//             if((i != Pstream::myProcNo()) &&  (procN == Pstream::myProcNo()))
//             {
//                 const int proc = i;
//                 {
//                     OPstream toNeighbour(proc);
//                     toNeighbour << cellsOnProcs[i].shrink() << molId 
//                                 << molPosition << cellI << cellId;
//                 }
//             }
//         }
//     }
// 
//     //- receiving
//     if(procN != Pstream::myProcNo())
//     {
//         List<label> cellLabels;
//         label idI;
//         vector posI;
//         label cellFromProc;
//         label cellIdFromProc;
// 
//         const int proc = procN;
//         {
//             IPstream fromNeighbour(proc);
//             fromNeighbour >> cellLabels >> idI >> posI 
//                           >> cellFromProc >> cellIdFromProc;
//         }
// 
//         vector rIJ;
//         scalar rIJMag;
//         scalar rIJMagSq;
//         scalar fIJMag;
//         label idJ;
//         scalar massJ;
// 
//         forAll(cellLabels, c)
//         {
//             const label& cell = cellLabels[c];
// 
//             forAll(cellOccupancy_[cell], mols)
//             {
//                 molecule* molJ = cellOccupancy_[cell][mols];
//         
//                 idJ = molJ->id();
//                 massJ = molJ->mass();
//                 rIJ = posI - molJ->position();
//                 
//                 rIJMagSq = magSqr(rIJ);
//                 
//                 if(pairPotentials_.rCutSqr(idI, idJ, rIJMagSq))
//                 {
//                     rIJMag = mag(rIJ);
//                 
//                     fIJMag = pairPotentials_.force(idI, idJ, rIJMag);
//                 
//                     // Acceleration increment for mol
//                     molJ->A() -= rIJ*fIJMag/(massJ);
//                 
//                     scalar potentialEnergy
//                     (
//                         pairPotentials_.energy(idI, idJ, rIJMag)
//                     );
//                 
//                     molJ->pE() -= 0.5*potentialEnergy;
//                 
//                     molJ->rDotf() -= 0.5*fIJMag*rIJMagSq;
//                 }
//             }
//         }
// 
//         //- rebuild referred cells
//         const labelList& cellRefIds = referredInteractionList_.refCellIds()[cellLabels[0]];
//         
//         forAll(cellRefIds, c)
//         {
//             waterReferredCell& refCellI = referredInteractionList_[cellRefIds[c]];
// 
//             if((refCellI.sourceProc() == procN) && (refCellI.sourceCell() == cellFromProc))
//             {
//                 DynamicList<waterReferredMolecule> referredMols(0);
// 
//                 forAll(refCellI, r)
//                 {
//                     if(r != cellIdFromProc)
//                     {
//                         referredMols.append
//                         (
//                             waterReferredMolecule(idI, posI)
//                         );
//                     }
//                 }
// 
//                 refCellI.clear();
//                 refCellI.transfer(referredMols.shrink());
//             }
//         }
//     }
// }


// void waterMoleculeCloud::insertMolInCellOccupancy(waterMolecule* mol)
// {
// 
// //     Info << "inserting molecule in cell: " << mol->cell() << endl;
// 
//     cellOccupancy_[mol->cell()].append(mol);
// 
//     cellOccupancy_[mol->cell()].shrink();
// }




// void waterMoleculeCloud::removeMolFromCellOccupancy
// (
//     waterMolecule* molI
// )
// {
//     DynamicList<waterMolecule*> updatedMolsInCell(0);
// 
//     const label& cellI = molI->cell();
// 
//     {
//         const List<waterMolecule*>& molsInCell = cellOccupancy_[cellI];
//     
//         forAll(molsInCell, m)
//         {
//             waterMolecule* molJ = molsInCell[m];
//     
//             if(molI != molJ)
//             {
//                 updatedMolsInCell.append(molJ);
//             }
//         }
//     }
// 
//     updatedMolsInCell.shrink();
//     cellOccupancy_[cellI].clear();
//     cellOccupancy_[cellI].transfer(updatedMolsInCell);
// }
// 
// 
// void waterMoleculeCloud::removeMolFromCellOccupancy
// (
//     const label& cellMolId,
//     const label& cell
// )
// {
//     DynamicList<waterMolecule*> molsInCell(0);
// 
//     forAll(cellOccupancy_[cell], c)
//     {
//         if(c != cellMolId)
//         {
//             molsInCell.append(cellOccupancy_[cell][c]);
//         }
//     }
// 
//     molsInCell.shrink();
//     cellOccupancy_[cell].clear();
//     cellOccupancy_[cell].transfer(molsInCell);
// }


label waterMoleculeCloud::numberOfMolsWithinRCut
(
    waterMolecule* mol
//     const potential& pot
)
{
    //- here we assume that the empty occupancy volume of the
    // neighbouring molecules is re-built in the next time-step.
    const pairPotentialList& pairPotentials_(pot_.pairPotentials());

    const waterReferredCellList& referredInteractionList_ = il_.ril();
    const waterDirectInteractionList& directInteractionCellList_ = il_.dil();
    const labelListList& fullInteractionCellList_ = directInteractionCellList_.fil();
//     const waterDirectInteractionList& fullInteractionCellList_ = il_.dil();

    const scalar& rCutMaxSqr_ = pairPotentials_.rCutMaxSqr();

    const label& cell = mol->cell();

//     const label& idI = mol->id();
    const vector& posI = mol->position();
//     const scalar& massI = mol->mass();

    vector rIJ;
    scalar rIJMagSq;
    //- real cells

    const labelList& dICL = fullInteractionCellList_[cell];

    label nMols = 0;

    //- molecules within direct interaction cells
    forAll(dICL, dCell)
    {
        List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            waterMolecule* molJ = molsCellJ[m];
            rIJ = posI - molJ->position();
            rIJMagSq = magSqr(rIJ);
            
            if(rIJMagSq <= rCutMaxSqr_)
            {
                nMols++;
            }
        }
    }

    //- molecules within cell
    forAll(cellOccupancy_[cell], mols)
    {
        waterMolecule* molJ = cellOccupancy_[cell][mols];

        if(molJ != mol)
        {
            rIJ = posI - molJ->position();
            rIJMagSq = magSqr(rIJ);
            
            if(rIJMagSq <= rCutMaxSqr_)
            {
                nMols++;
            }
        }
    }

    const labelList& refCellIds = referredInteractionList_.refCellIds()[cell];

    forAll(refCellIds, r)
    {
        const label& refCellId = refCellIds[r];
        const waterReferredCell& refCellI = referredInteractionList_[refCellId];

        forAll(refCellI, refMols)
        {
            const waterReferredMolecule* molJ = &(refCellI[refMols]);

            rIJ = posI - molJ->position();
            rIJMagSq = magSqr(rIJ);
            
            if(rIJMagSq <= rCutMaxSqr_)
            {
                nMols++;
            }
        }
    }

    return nMols;

}


void waterMoleculeCloud::updateNeighbouringRadii
(
    waterMolecule* molI
)
{
    molI->R() = GREAT;
    const label& cell = molI->cell();
//     const label& idI = mol->id();
//     const waterMolecule::constantProperties& constPropI = constProps(idI);

//     const pairPotentialList& pairPot(pot_.pairPotentials());

    const labelListList& fullInteractionCellList = il_.dil().fil();
//     const pairPotential& electrostatic = pairPot.electrostatic();

    //- real cells

//     label idJ = -1;

    const labelList& dICL = fullInteractionCellList[cell];

    //- molecules within direct interaction cells (not incl. owner cell)
    forAll(dICL, dCell)
    {
        const List< waterMolecule* > molsCellJ = cellOccupancy_[dICL[dCell]];

        forAll(molsCellJ, m)
        {
            waterMolecule* molJ = molsCellJ[m];

            vector rIJ = molI->position() - molJ->position();
    
            scalar rIJMag = mag(rIJ);
    
            if(molI->R() > rIJMag)
            {
                molI->R() = rIJMag;
            }
    
            if(molJ->R() > rIJMag)
            {
                molJ->R() = rIJMag;
            }
        }
    }

    //- molecules within owner cell

    forAll(cellOccupancy_[cell], mols)
    {
        waterMolecule* molJ = cellOccupancy_[cell][mols];

        if(molJ != molI)
        {
            vector rIJ = molI->position() - molJ->position();
    
            scalar rIJMag = mag(rIJ);
    
            if(molI->R() > rIJMag)
            {
                molI->R() = rIJMag;
            }
    
            if(molJ->R() > rIJMag)
            {
                molJ->R() = rIJMag;
            }
        }
    }

    //- referred cells 

    // refer in molecules due to periodic boundary conditions
    {
        waterReferredCellList& ril = il_.ril();
    
        const labelList& refCellsPerIds = ril.refCellIdsPer()[molI->cell()];
    
        forAll(refCellsPerIds, r)
        {
            const label& refCellId = refCellsPerIds[r];
    
            waterReferredCell& refCellI = ril[refCellId];//*****
    
            refCellI.referInMolecule
            (
                waterReferredMolecule
                (
                    molI->id(),
                    molI->position(),
                    molI->sitePositions(),
                    molI->v(),
                    molI->fraction(),
                    molI->trackingNumber()
                )
            );
        }
    }

    const waterReferredCellList& referredInteractionList = il_.ril();
    const labelList& refCellIds = referredInteractionList.refCellIds()[cell];

    forAll(refCellIds, r)
    {
        const label& refCellId = refCellIds[r];
        const waterReferredCell& refCellI = referredInteractionList[refCellId];

        forAll(refCellI, refMols)
        {
            const waterReferredMolecule* molRef = &(refCellI[refMols]);

            vector rRealRef = molI->position() - molRef->position();
    
            scalar rRealRefMag = mag(rRealRef);
    
            if(molI->R() > rRealRefMag)
            {
                molI->R() = rRealRefMag;
            }
        }
    }
}

void waterMoleculeCloud::updateRadii()
{
//     const pairPotentialList& pairPot(pot_.pairPotentials());
//     const pairPotential& electrostatic = pairPot.electrostatic();

    const waterReferredCellList& referredInteractionList = il_.ril();
    
    forAll(referredInteractionList, refCellId)
    {
        const waterReferredCell& refCellI = referredInteractionList[refCellId];

        forAll(refCellI, refMols)
        {
            const waterReferredMolecule* molRef = &(refCellI[refMols]);

//             label idRef = molRef->id();

       
            const labelList& realCells = refCellI.realCellsForInteraction();
            
            forAll(realCells, rC)
            {
                const label& realCellI = realCells[rC];
    
                List<waterMolecule*>& molsInCell = cellOccupancy_[realCellI];
    
                forAll(molsInCell, m)
                {
                    waterMolecule* molReal = molsInCell[m];

                    vector rRealRef = molReal->position() - molRef->position();
            
                    scalar rRealRefMag = mag(rRealRef);
            
                    if(molReal->R() > rRealRefMag)
                    {
                        molReal->R() = rRealRefMag;
                    }
                }
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
