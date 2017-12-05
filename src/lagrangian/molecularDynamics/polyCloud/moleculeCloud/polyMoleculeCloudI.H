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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// NEW
inline void Foam::polyMoleculeCloud::evaluatePair
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
  
    label idI = molI->id();
    label idJ = molJ->id();

    // pair potential interactions
    controllers_.controlDuringForceComputation(molI, molJ); 

    fields_.measurementsDuringForceComputation
    (
        molI,
        molJ
    );

    if(!molI->frozen() || !molJ->frozen())
    {
        // fraction
        scalar f = molI->fraction();
    
        if(molJ->fraction() < f)
        {
            f = molJ->fraction();
        }
        
        // looping over just pair potential sites
        forAll(cP_.pairSitesLinks()[idI], i)
        {
            // sI is the index of the full list of sites
            label sI = cP_.pairSitesLinks()[idI][i];
            
            forAll(cP_.pairSitesLinks()[idJ], j)
            {
                // sJ is the index of the full list of sites
                label sJ = cP_.pairSitesLinks()[idJ][j];

                vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);
                
                // determine index of pairPotential from list of pairPotentials_
                
//                 label k = p_.pairPots().getIndexFromPairPotentialSites(idI, idJ, i, j);
                
                label k = p_.pairPots().pairPotentialIndex(idI, idJ, sI, sJ);
                
//                 Info << nl << "idI = " << idI << " = " << cP_.siteNames(idI) 
//                      << ", idJ = " << idJ << " = " << cP_.siteNames(idJ)
//                      << endl;
//                      
//                 Info << "posI = " << molI->position() 
//                      << ", posJ = " << molJ->position() 
//                      << endl;
//                      
//                 Info << "sIdI = " << sI << ", sIdJ = " << sJ 
//                      << ", siteI = " << molI->sitePositions()[sI]
//                      << ", siteJ = " << molJ->sitePositions()[sJ]
//                      << ", mag = " << mag(rsIsJ)        
//                      << endl;                
//                      
//                 Info << "k = " << k << endl;
                     
                if(p_.pairPots().rCutSqr(k, rsIsJMagSq))
                {
//                     if(!p_.pairPots().excludeSites(molI, molJ, sI, sJ, k))
                    {
                        scalar rsIsJMag = mag(rsIsJ);

                        vector fsIsJ = f*(rsIsJ/rsIsJMag)*p_.pairPots().force(k, rsIsJMag);

//                         Info << "force mag = "<< p_.pairPots().force(k, rsIsJMag) << endl;
//                         Info << "force = " << fsIsJ << endl;
                             
                        molI->siteForces()[sI] += fsIsJ;
                        molJ->siteForces()[sJ] += -fsIsJ;

                        scalar potentialEnergy
                        (
                            f*p_.pairPots().energy(k, rsIsJMag)
                        );

                        molI->potentialEnergy() += 0.5*potentialEnergy;
                        molJ->potentialEnergy() += 0.5*potentialEnergy;
            
                        vector rIJ = molI->position() - molJ->position();
            
                        tensor virialContribution = (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
            
                        molI->rf() += virialContribution;
                        molJ->rf() += virialContribution;

                        fields_.measurementsDuringForceComputationSite
                        (
                            molI,
                            molJ,
                            sI,
                            sJ
                        );
                    }
                }
            }
        }
 
        {
            Foam::vector rIJ = molI->position() - molJ->position();
    
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

        
        //chargePotNamesToSites
        forAll(cP_.chargeSitesLinks()[idI], pI)
        {
            label sI = cP_.chargeSitesLinks()[idI][pI];
    
            forAll(cP_.chargeSitesLinks()[idJ], pJ)
            {
                label sJ = cP_.chargeSitesLinks()[idJ][pJ];
    
                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                scalar rsIsJMagSq = magSqr(rsIsJ);

/*                Info << "idI = " << idI << " = " << cP_.siteNames(idI) 
                     << ", idJ = " << idJ << " = " << cP_.siteNames(idJ)
                     << endl;
                     
                Info << "posI = " << molI->position() 
                     << ", posJ = " << molJ->position() 
                     << endl;
                     
                Info << "sIdI = " << sI << ", sIdJ = " << sJ 
                     << ", siteI = " << molI->sitePositions()[sI]
                     << ", siteJ = " << molJ->sitePositions()[sJ]
                     << ", mag = " << mag(rsIsJ)        
                     << endl;    */                   
                //if(rsIsJMagSq <= p_.pairPots().electrostaticPot()->rCutSqr())
                if(rsIsJMagSq < p_.pairPots().electrostaticPot()->rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    scalar chargeI = cP_.siteCharges(idI)[sI];
                    scalar chargeJ = cP_.siteCharges(idJ)[sJ];
        
                    vector fsIsJ =
                        f*(rsIsJ/rsIsJMag)
                        *chargeI*chargeJ
                        *p_.pairPots().electrostaticPot()->force(rsIsJMag);

//                     Info << "force = " << fsIsJ << endl;                        
                        
                    molI->siteForces()[sI] += fsIsJ;
        
                    molJ->siteForces()[sJ] += -fsIsJ;
        
                    scalar potentialEnergy =
                        f*chargeI*chargeJ
                        *p_.pairPots().electrostaticPot()->energy(rsIsJMag);
        
                    molI->potentialEnergy() += 0.5*potentialEnergy;
        
                    molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                    vector rIJ = molI->position() - molJ->position();
        
                    tensor virialContribution =
                        (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
        
                    molI->rf() += virialContribution;
                    molJ->rf() += virialContribution;
    
                    fields_.measurementsDuringForceComputationSite
                    (
                        molI,
                        molJ,
                        sI,
                        sJ
                    );
                }
            }
        }
    }
}


inline bool Foam::polyMoleculeCloud::evaluatePotentialLimit
(
    polyMolecule* molI,
    polyMolecule* molJ,
    const scalar& potentialEnergyLimit
) 
{
    label idI = molI->id();
    label idJ = molJ->id();

    if(!molI->frozen() || !molJ->frozen())
    {
        // fraction
        scalar f = molI->fraction();
    
        if(molJ->fraction() < f)
        {
            f = molJ->fraction();
        }
        
        // looping over just pair potential sites
        forAll(cP_.pairSitesLinks()[idI], i)
        {
            // sI is the index of the full list of sites
            label sI = cP_.pairSitesLinks()[idI][i];
            
            forAll(cP_.pairSitesLinks()[idJ], j)
            {
                // sJ is the index of the full list of sites
                label sJ = cP_.pairSitesLinks()[idJ][j];

                vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);
                
                label k = p_.pairPots().pairPotentialIndex(idI, idJ, sI, sJ);
                
                if(p_.pairPots().rCutSqr(k, rsIsJMagSq))
                {
                    if(!p_.pairPots().excludeSites(molI, molJ, sI, sJ, k))
                    {
                        scalar rsIsJMag = mag(rsIsJ);

                        scalar rMin = p_.pairPots().rMin(k);
                        
                        // Guard against pairPot.energy being evaluated
                        // if rIJMag < SMALL. A floating point exception will
                        // happen otherwise.
                        
                        
                        if (rsIsJMag < SMALL)
                        {
                            return true;
                        }
            
                        // Guard against pairPot.energy being evaluated if rIJMag <
                        // rMin.  A tabulation lookup error will occur otherwise.
            
                        if (rsIsJMag < rMin)
                        {
                            return true;
                        }
        
                        if
                        (
                            mag(p_.pairPots().energy(k, rsIsJMag))
                            > potentialEnergyLimit
                        )
                        {
                            return true;
                        }
                    }
                }
            }
        }
        
        forAll(cP_.chargeSitesLinks()[idI], pI)
        {
            label sI = cP_.chargeSitesLinks()[idI][pI];
    
            forAll(cP_.chargeSitesLinks()[idJ], pJ)
            {
                label sJ = cP_.chargeSitesLinks()[idJ][pJ];
    
                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                scalar rsIsJMagSq = magSqr(rsIsJ);

                if(rsIsJMagSq < p_.pairPots().electrostaticPot()->rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);
                    
                    // Guard against pairPot.energy being evaluated
                    // if rIJMag < SMALL. A floating point exception will
                    // happen otherwise.
    
                    if (rsIsJMag < SMALL)
                    {
                        return true;
                    }
    
                    if (rsIsJMag <  p_.pairPots().electrostaticPot()->rMin())
                    {
                        return true;
                    }
                    
                    scalar chargeI = cP_.siteCharges(idI)[sI];
                    scalar chargeJ = cP_.siteCharges(idJ)[sJ];
    
                    if
                    (
                        mag(chargeI*chargeJ*p_.pairPots().electrostaticPot()->energy(rsIsJMag))
                        > potentialEnergyLimit
                    )
                    {
                        return true;
                    }
                }
            }
        }
    }
    
    return false;
}

/*
// OLD
inline void Foam::polyMoleculeCloud::evaluatePair
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
    const pairPotentialList& pairPot = pot_.pairPotentials();
    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI->id();
    label idJ = molJ->id();

    const polyMolecule::constantProperties& constPropI(constProps(idI));
    const polyMolecule::constantProperties& constPropJ(constProps(idJ));

    // pair potential interactions
    controllers_.controlDuringForceComputation(molI, molJ); 

    fields_.measurementsDuringForceComputation
    (
        molI,
        molJ
    );

    if(!molI->frozen() || !molJ->frozen())
    {
        // fraction
        scalar f = molI->fraction();
    
        if(molJ->fraction() < f)
        {
            f = molJ->fraction();
        }

        forAll(constPropI.pairPotSites(), pI)
        {
        	label sI = constPropI.pairPotSites()[pI];

            forAll(constPropJ.pairPotSites(), pJ)
            {
                label sJ = constPropJ.pairPotSites()[pJ];

                Foam::vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                label idsI = constPropI.sites()[sI].siteId();
                label idsJ = constPropJ.sites()[sJ].siteId();

                if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    Foam::vector fsIsJ = (rsIsJMag == 0.0? Foam::vector::zero : (f * (rsIsJ/rsIsJMag) * pairPot.force(idsI, idsJ, rsIsJMag)));

                    molI->siteForces()[sI] += fsIsJ;

                    molJ->siteForces()[sJ] += -fsIsJ;

//                     Info<< "position I = " << molI->position() 
//                         << ", site = " << idsI 
//                         << ", mag = " << rsIsJMag
//                         << ", force = " << fsIsJ 
//                         << endl;
                        
                    scalar potentialEnergy
                    (
                        f*pairPot.energy(idsI, idsJ, rsIsJMag)
                    );

                    molI->potentialEnergy() += 0.5*potentialEnergy;
        
                    molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                    Foam::vector rIJ = molI->position() - molJ->position();
        
                    Foam::tensor virialContribution = (rsIsJMagSq == 0.0? (rsIsJ*fsIsJ)*(rsIsJ & rIJ): (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq);
        
                    molI->rf() += virialContribution;
                    molJ->rf() += virialContribution;

                    fields_.measurementsDuringForceComputationSite
                    (
                        molI,
                        molJ,
                        sI,
                        sJ
                    );
                }
            }
        }

        {
            Foam::vector rIJ = molI->position() - molJ->position();
    
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

        forAll(constPropI.electrostaticSites(), pI)
        {
            label sI = constPropI.electrostaticSites()[pI];
    
            forAll(constPropJ.electrostaticSites(), pJ)
            {
                label sJ = constPropJ.electrostaticSites()[pJ];
    
                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                scalar rsIsJMagSq = magSqr(rsIsJ);
        
                if(rsIsJMagSq < electrostatic.rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);
        
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
        
                    vector fsIsJ =
                        f*(rsIsJ/rsIsJMag)
                        *chargeI*chargeJ*electrostatic.force(rsIsJMag);
        
                    molI->siteForces()[sI] += fsIsJ;
        
                    molJ->siteForces()[sJ] += -fsIsJ;
        
                    scalar potentialEnergy =
                        f*chargeI*chargeJ
                        *electrostatic.energy(rsIsJMag);
        
                    molI->potentialEnergy() += 0.5*potentialEnergy;
        
                    molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                    vector rIJ = molI->position() - molJ->position();
        
                    tensor virialContribution =
                        (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
        
                    molI->rf() += virialContribution;
                    molJ->rf() += virialContribution;
    
                    fields_.measurementsDuringForceComputationSite
                    (
                        molI,
                        molJ,
                        sI,
                        sJ
                    );
                }
            }
        }
    }
}*/

/*
inline bool Foam::polyMoleculeCloud::evaluatePotentialLimit
(
    polyMolecule* molI,
    polyMolecule* molJ,
    const scalar& potentialEnergyLimit
) const
{
    const pairPotentialList& pairPot = pot_.pairPotentials();
    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI->id();
    label idJ = molJ->id();

    const polyMolecule::constantProperties& constPropI(constProps(idI));
    const polyMolecule::constantProperties& constPropJ(constProps(idJ));

    if(!molI->frozen() || !molJ->frozen())
    {
        forAll(constPropI.pairPotSites(), pI)
        {
            label sI = constPropI.pairPotSites()[pI];
    
            forAll(constPropJ.pairPotSites(), pJ)
            {
                label sJ = constPropJ.pairPotSites()[pJ];
    
                label idsI = constPropI.sites()[sI].siteId();
                label idsJ = constPropJ.sites()[sJ].siteId();
        
                vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
        
                scalar rsIsJMagSq = magSqr(rsIsJ);
        
                if (pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);
        
                    // Guard against pairPot.energy being evaluated
                    // if rIJMag < SMALL. A floating point exception will
                    // happen otherwise.
        
                    if (rsIsJMag < SMALL)
                    {
                        return true;
                    }
        
                    // Guard against pairPot.energy being evaluated if rIJMag <
                    // rMin.  A tabulation lookup error will occur otherwise.
        
                    if (rsIsJMag < pairPot.rMin(idsI, idsJ))
                    {
                        return true;
                    }
    
                    if
                    (
                        mag(pairPot.energy(idsI, idsJ, rsIsJMag))
                        > potentialEnergyLimit
                    )
                    {
                        return true;
                    };
                }
            }
        }
    
        forAll(constPropI.electrostaticSites(), pI)
        {
            label sI = constPropI.electrostaticSites()[pI];
    
            forAll(constPropJ.electrostaticSites(), pJ)
            {
                label sJ = constPropJ.electrostaticSites()[pJ];
    
                vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                scalar rsIsJMagSq = magSqr(rsIsJ);
    
                
                if(rsIsJMagSq < electrostatic.rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);
    
                    // Guard against pairPot.energy being evaluated
                    // if rIJMag < SMALL. A floating point exception will
                    // happen otherwise.
    
                    if (rsIsJMag < SMALL)
                    {
                        return true;
                    }
    
                    if (rsIsJMag < electrostatic.rMin())
                    {
                        return true;
                    }
    
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
    
                    if
                    (
                        mag(chargeI*chargeJ*electrostatic.energy(rsIsJMag))
                        > potentialEnergyLimit
                    )
                    {
                        return true;
                    };
                }
            }
        }
    }
    return false;
}
*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::polyMoleculeCloud::mesh() const
{
    return mesh_;
}

// inline const Foam::potential& Foam::polyMoleculeCloud::pot() const
// {
//     return pot_;
// }

inline const Foam::List<Foam::DynamicList<Foam::polyMolecule*> >&
    Foam::polyMoleculeCloud::cellOccupancy() const
{
    return cellOccupancy_;
}

inline const Foam::potentials& Foam::polyMoleculeCloud::pot() const
{
    return p_;
}

inline Foam::potentials& Foam::polyMoleculeCloud::pot()
{
    return p_;
}


inline const Foam::constantMoleculeProperties& Foam::polyMoleculeCloud::cP() const
{
    return cP_;
}

// inline const Foam::List<Foam::polyMolecule::constantProperties>
//     Foam::polyMoleculeCloud::constProps() const
// {
//     return constPropList_;
// }
// 
// inline const Foam::polyMolecule::constantProperties&
//     Foam::polyMoleculeCloud::constProps(label id) const
// {
//     return constPropList_[id];
// }

inline Foam::cachedRandomMD& Foam::polyMoleculeCloud::rndGen()
{
    return rndGen_;
}

inline const Foam::reducedUnits& Foam::polyMoleculeCloud::redUnits() const
{
    return redUnits_;
}

inline const Foam::polyFieldProperties& Foam::polyMoleculeCloud::fields() const
{
    return fields_;
}

inline const  Foam::cyclicBoundaries& Foam::polyMoleculeCloud::cyclics() const
{
    return cyclics_;
}

inline Foam::polyFieldProperties& Foam::polyMoleculeCloud::fields()
{
    return fields_;
}

inline const Foam::polyBoundaries& Foam::polyMoleculeCloud::boundaries() const
{
    return boundaries_;
}

inline Foam::polyBoundaries& Foam::polyMoleculeCloud::boundaries()
{
    return boundaries_;
}

inline const Foam::polyControllers& Foam::polyMoleculeCloud::controllersInfo() const
{
    return controllers_;
}

inline Foam::polyControllers& Foam::polyMoleculeCloud::controllersInfo()
{
    return controllers_;
}

inline const Foam::polyFaceTracker& Foam::polyMoleculeCloud::tracker() const
{
    return trackingInfo_;
}

inline Foam::polyFaceTracker& Foam::polyMoleculeCloud::tracker()
{
    return trackingInfo_;
}

inline Foam::cellInteractions<Foam::polyMolecule>& Foam::polyMoleculeCloud::il()
{
    return iL_;
}

inline Foam::trackingNumber& Foam::polyMoleculeCloud::moleculeTracking()
{
    return moleculeTracking_;
}

inline Foam::clockTimer& Foam::polyMoleculeCloud::clock()
{
    return clock_;
}

// ************************************************************************* //
