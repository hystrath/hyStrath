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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
inline Foam::polyMolecule::constantProperties::constantProperties()
:
    nSites_(0),
    sites_(),
    electrostaticSites_(),
    pairPotSites_(),
    momentOfInertia_(diagTensor(0, 0, 0)),
    mass_(0.0)
{}


inline Foam::polyMolecule::constantProperties::constantProperties
(
    const dictionary& dict,
    const reducedUnits& rU,
    const labelList& siteIds
)
:
    nSites_(siteIds.size()),
    sites_(siteIds.size()),
    electrostaticSites_(),
    pairPotSites_(),
    momentOfInertia_(),
    mass_(0.0)
{
    // read properties
    Field<vector> siteReferencePositions(dict.lookup("siteReferencePositions"));
    List<scalar> siteMasses(dict.lookup("siteMasses"));
    List<scalar> siteCharges(dict.lookup("siteCharges"));
    List<word> idNames(dict.lookup("siteIds"));
    List<word> pairPotentialIds(dict.lookup("pairPotentialSiteIds"));

    // check

    if
    (
        ( siteMasses.size() != nSites_ )
     || ( idNames.size() != nSites_ )
     || ( siteReferencePositions.size() != nSites_ )
     || ( siteCharges.size() != nSites_ )
    )
    {
        FatalErrorIn("polyMolecule::constantProperties::constantProperties")
            << "Sizes of site masses, charges and "
            << "referencePositions are not the same as ids: " << nSites_
            << nl << abort(FatalError);
    }

    // set each site

    label nChargedMols = 0;
    label nPairPotMols = 0;

    forAll(sites_, s)
    {
        bool charge = false;
        bool pairPot = false;

        if(mag(siteCharges[s]) > VSMALL)
        {
            charge = true;
            nChargedMols++;
        }

        if( findIndex(pairPotentialIds, idNames[s]) != -1 )
        {
            pairPot = true;
            nPairPotMols++;
        }

        if(!charge && !pairPot)
        {
            Info<< "WARNING: " << idNames[s] << " is a non-interacting site."
                << endl;
        }

        sites_[s] = constPropSite
                    (
                        siteReferencePositions[s],
                        siteMasses[s],
                        siteCharges[s],
                        siteIds[s],
                        idNames[s],
                        pairPot,
                        charge
                    );
    }

    // set the electrostatic and pairPotential site lists
    electrostaticSites_.setSize(nChargedMols, -1);
    pairPotSites_.setSize(nPairPotMols, -1);

    label counterI = 0;
    label counterJ = 0;

    forAll(sites_, s)
    {
        if(sites_[s].pairPotentialSite())
        {
            pairPotSites_[counterI] = s;
            counterI++; 
        }

        if(sites_[s].electrostaticSite())
        {
            electrostaticSites_[counterJ] = s;
            counterJ++; 
        }
    }

    Info << "site ids: " << siteIds << endl;
    Info << "pairPotSites_: " << pairPotSites_ << endl;
    Info << "electrostaticSites_: " << electrostaticSites_ << endl;

    // reduced units

    if(rU.runReducedUnits())
    {
        forAll(sites_, s)
        {
            sites_[s].siteCharge() /= rU.refCharge();
            sites_[s].siteMass() /= rU.refMass();
            sites_[s].siteReferencePosition() /= rU.refLength();
            mass_ += sites_[s].siteMass();
        }
    }

    vector centreOfMass(vector::zero);

    // Calculate the centre of mass of the body and subtract it from each
    // position

    forAll(sites_, s)
    {
        centreOfMass += sites_[s].siteReferencePosition()*sites_[s].siteMass();
    }

    centreOfMass /= mass_;

    forAll(sites_, s)
    {
        sites_[s].siteReferencePosition() -= centreOfMass;
    }

    if (sites_.size() == 1)
    {
        sites_[0].siteReferencePosition() = vector::zero;

        momentOfInertia_ = diagTensor(-1, -1, -1);

//         FatalErrorIn("polyMolecule::constantProperties::constantProperties")
//             << "Use the atomistic cloud "
//             << ", since no of sites = " << sites_.size()
//             << nl << abort(FatalError);

    }
    else if (linearMoleculeTest())
    {
        // Linear molecule.

        Info<< nl << "Linear molecule." << endl;

        vector dir = sites_[1].siteReferencePosition() - sites_[0].siteReferencePosition();

        dir /= mag(dir);

        tensor Q = rotationTensor(dir, vector(1,0,0));

        // Transform the site positions
        forAll(sites_, s)
        {
            sites_[s].siteReferencePosition() = (Q & sites_[s].siteReferencePosition());
        }


        // The rotation was around the centre of mass but remove any
        // components that have crept in due to floating point errors

        centreOfMass = vector::zero;

        forAll(sites_, s)
        {
            centreOfMass += sites_[s].siteReferencePosition()*sites_[s].siteMass();
        }

        centreOfMass /= mass_;

        forAll(sites_, s)
        {
            sites_[s].siteReferencePosition() -= centreOfMass;
        }

        diagTensor momOfInertia = diagTensor::zero;

        forAll(sites_, s)
        {
            const vector& p(sites_[s].siteReferencePosition());
    
            momOfInertia += sites_[s].siteMass()
                *diagTensor(0, p.x()*p.x(), p.x()*p.x());
        }

        momentOfInertia_ = diagTensor
        (
            -1,
            momOfInertia.yy(),
            momOfInertia.zz()
        );
    }
    else
    {
        // Fully 6DOF molecule

        // Calculate the inertia tensor in the current orientation

        tensor momOfInertia(tensor::zero);
    
        forAll(sites_, s)
        {
            const vector& p(sites_[s].siteReferencePosition());
    
            momOfInertia += sites_[s].siteMass()*tensor
            (
                p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
            );
        }
        
//         Info << "Eigen values: x = " << eigenValues(momOfInertia).x()
//             << " y = " << eigenValues(momOfInertia).y()
//             << " z = " << eigenValues(momOfInertia).z()
//             << endl;
        Info << "momOfInertia (full tensor) = " <<  momOfInertia << endl;
        
        if (eigenValues(momOfInertia).x() < VSMALL)
        {
//             FatalErrorIn("polyMolecule::constantProperties::constantProperties")
//                 << "An eigenvalue of the inertia tensor is zero, the molecule "
//                 << dict.name() << " is not a valid 6DOF shape."
//                 << nl << abort(FatalError);
            
            Info << "Warning no adjustment to be made to molecule" << endl;
            
            momentOfInertia_ = diagTensor
            (
                momOfInertia.xx(),
                momOfInertia.yy(),
                momOfInertia.zz()
            ); 
        }
        else
        {
            Info << "Adjusting molecule" << endl;
            
            // Normalise the inertia tensor magnitude to avoid SMALL numbers in the
            // components causing problems
    
            momOfInertia /= eigenValues(momOfInertia).x();
        
            tensor e = eigenVectors(momOfInertia);
    
        
            // Calculate the transformation between the principle axes to the
            // global axes
        
            tensor Q = vector(1,0,0)*e.x() + vector(0,1,0)*e.y() + vector(0,0,1)*e.z();
        
        //         Info << "Q: " << Q << endl;
        
            // Transform the site positions
            forAll(sites_, s)
            {
                sites_[s].siteReferencePosition() = (Q & sites_[s].siteReferencePosition());
            }

            // Recalculating the moment of inertia with the new site positions
        
            // The rotation was around the centre of mass but remove any
            // components that have crept in due to floating point errors
        
            centreOfMass = vector::zero;
        
            forAll(sites_, s)
            {
                centreOfMass += sites_[s].siteReferencePosition()*sites_[s].siteMass();
            }
        
            centreOfMass /= mass_;
        
            forAll(sites_, s)
            {
                sites_[s].siteReferencePosition() -= centreOfMass;

            }
            
            // Calculate the moment of inertia in the principle component
            // reference frame
        
            momOfInertia = tensor::zero;
        
            forAll(sites_, s)
            {
                const vector& p(sites_[s].siteReferencePosition());
        
                momOfInertia += sites_[s].siteMass()*tensor
                (
                    p.y()*p.y() + p.z()*p.z(), -p.x()*p.y(), -p.x()*p.z(),
                    -p.y()*p.x(), p.x()*p.x() + p.z()*p.z(), -p.y()*p.z(),
                    -p.z()*p.x(), -p.z()*p.y(), p.x()*p.x() + p.y()*p.y()
                );
            }
        
            momentOfInertia_ = diagTensor
            (
                momOfInertia.xx(),
                momOfInertia.yy(),
                momOfInertia.zz()
            );            
        }
        
        Info << "moment of inertia (diag tensor): " << momentOfInertia_ << endl;
    }
}*/


inline Foam::polyMolecule::polyMolecule
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    const vector& specialPosition,
//     const constantProperties& constProps,
    const constantMoleculeProperties& cP,
    const label special,
    const label id,
    const scalar& frac,
    const label trackingNumber
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    Q_(Q),
    v_(v),
    a_(a),
    pi_(pi),
    tau_(tau),
    specialPosition_(specialPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(special),
    id_(id),
    R_(GREAT),
    frac_(frac),
    trackingNumber_(trackingNumber),
//     siteForces_(constProps.nSites(), vector::zero),
//     sitePositions_(constProps.nSites())
    siteForces_(cP.siteNames(id).size(), vector::zero),
    sitePositions_(cP.siteNames(id).size())
{
    setSitePositions(cP);
}

/*
inline bool Foam::polyMolecule::constantProperties::linearMoleculeTest() const
{
    if (sites_.size() == 2)
    {
        return true;
    }

    vector refDir = sites_[1].siteReferencePosition() - sites_[0].siteReferencePosition();

    refDir /= mag(refDir);

    for
    (
        label i = 2;
        i < nSites_;
        i++
    )
    {
        vector dir = sites_[i].siteReferencePosition() - sites_[i-1].siteReferencePosition();

        dir /= mag(dir);

        if (mag(refDir & dir) < 1 - SMALL)
        {
            return false;
        }
    }

    return true;
}*/


// * * * * * * * constantProperties Member Functions * * * * * * * * * * * * //
// inline const Foam::List< Foam::constPropSite>&
// Foam::polyMolecule::constantProperties::sites() const
// {
//     return sites_;
// }

/*
inline const Foam::List<Foam::label>&
Foam::polyMolecule::constantProperties::pairPotSites() const
{
    return pairPotSites_;
}

inline const Foam::List<Foam::label>&
Foam::polyMolecule::constantProperties::electrostaticSites() const
{
    return electrostaticSites_;
}

inline const Foam::diagTensor&
Foam::polyMolecule::constantProperties::momentOfInertia() const
{
    return momentOfInertia_;
}


inline bool Foam::polyMolecule::constantProperties::linearMolecule() const
{
    return ((momentOfInertia_.xx() < 0) && (momentOfInertia_.yy() > 0));
}


inline bool Foam::polyMolecule::constantProperties::pointMolecule() const
{
    return (momentOfInertia_.zz() < 0);
}


inline Foam::label Foam::polyMolecule::constantProperties::degreesOfFreedom() const
{
    if (linearMolecule())
    {
        return 5;
    }
    else if (pointMolecule())
    {
        return 3;
    }
    else
    {
        return 6;
    }
}*/

/*
inline Foam::scalar Foam::polyMolecule::constantProperties::mass() const
{
    return mass_;
}


inline Foam::label Foam::polyMolecule::constantProperties::nSites() const
{
    return sites_.size();

}*/

// * * * * * * * * * * * * polyMolecule Member Functions * * * * * * * * * * * * //

inline const Foam::tensor& Foam::polyMolecule::Q() const
{
    return Q_;
}


inline Foam::tensor& Foam::polyMolecule::Q()
{
    return Q_;
}


inline const Foam::vector& Foam::polyMolecule::v() const
{
    return v_;
}


inline Foam::vector& Foam::polyMolecule::v()
{
    return v_;
}


inline const Foam::vector& Foam::polyMolecule::a() const
{
    return a_;
}


inline Foam::vector& Foam::polyMolecule::a()
{
    return a_;
}


inline const Foam::vector& Foam::polyMolecule::pi() const
{
    return pi_;
}


inline Foam::vector& Foam::polyMolecule::pi()
{
    return pi_;
}


inline const Foam::vector& Foam::polyMolecule::tau() const
{
    return tau_;
}


inline Foam::vector& Foam::polyMolecule::tau()
{
    return tau_;
}


inline const Foam::List<Foam::vector>& Foam::polyMolecule::siteForces() const
{
    return siteForces_;
}


inline Foam::List<Foam::vector>& Foam::polyMolecule::siteForces()
{
    return siteForces_;
}


inline const Foam::List<Foam::vector>& Foam::polyMolecule::sitePositions() const
{
    return sitePositions_;
}


inline Foam::List<Foam::vector>& Foam::polyMolecule::sitePositions()
{
    return sitePositions_;
}


inline const Foam::vector& Foam::polyMolecule::specialPosition() const
{
    return specialPosition_;
}


inline Foam::vector& Foam::polyMolecule::specialPosition()
{
    return specialPosition_;
}


inline Foam::scalar Foam::polyMolecule::potentialEnergy() const
{
    return potentialEnergy_;
}


inline Foam::scalar& Foam::polyMolecule::potentialEnergy()
{
    return potentialEnergy_;
}


inline const Foam::tensor& Foam::polyMolecule::rf() const
{
    return rf_;
}


inline Foam::tensor& Foam::polyMolecule::rf()
{
    return rf_;
}

inline Foam::label Foam::polyMolecule::special() const
{
    return special_;
}

inline Foam::label& Foam::polyMolecule::special()
{
    return special_;
}

inline bool Foam::polyMolecule::tethered() const
{
    return special_ == SPECIAL_TETHERED;
}

inline Foam::label Foam::polyMolecule::id() const
{
    return id_;
}

inline Foam::label& Foam::polyMolecule::id()
{
    return id_;
}

inline const Foam::scalar& Foam::polyMolecule::R() const
{
    return R_;
}

inline Foam::scalar& Foam::polyMolecule::R()
{
    return R_;
}

inline const Foam::scalar& Foam::polyMolecule::fraction() const
{
    return frac_;
}

inline Foam::scalar& Foam::polyMolecule::fraction()
{
    return frac_;
}

inline const Foam::label& Foam::polyMolecule::trackingNumber() const
{
    return trackingNumber_;
}

inline Foam::label& Foam::polyMolecule::trackingNumber()
{
    return trackingNumber_;
}

inline bool Foam::polyMolecule::frozen()
{
    if(special_ == SPECIAL_FROZEN)
    {
        return true;
    }
    else
    {
        return false;
    }
}

inline bool Foam::polyMolecule::referred()
{
    return special_ == SPECIAL_REFERRED;
}


// ************************************************************************* //
