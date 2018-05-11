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
    dsmcSpherical

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcSpherical.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcSpherical, 0);

    addToRunTimeSelectionTable
    (
        dsmcCoordinateSystem, 
        dsmcSpherical,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcSpherical::sphericalWeighting()
{
    forAll(cloud_.cellOccupancy(), c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cloud_.cellOccupancy()[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            
            const scalar oldRadialWeight = p->RWF();
                        
            const scalar newRadialWeight = RWF(c);

            p->RWF() = newRadialWeight;
            
            if (oldRadialWeight > newRadialWeight) 
            {
                //- particle might be cloned
                scalar prob = (oldRadialWeight/newRadialWeight) - 1.0;
                
                while(prob > 1.0)
                {
                    //- add a particle and reduce prob by 1.0
                    vector U = p->U();
                    
                    cloud_.addNewParcel
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
                
                if (prob > cloud_.rndGen().scalar01())
                {
                    vector U = p->U();
                    
                    cloud_.addNewParcel
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
            else if (newRadialWeight > oldRadialWeight)
            {           
                //- particle might be deleted
                if ((oldRadialWeight/newRadialWeight) < cloud_.rndGen().scalar01())
                {
                    cloud_.deleteParticle(*p);
                } 
            } 
        }
    }
}


void dsmcSpherical::updateRWF()
{
    forAll(RWF_, celli)
    {
        RWF_[celli] = recalculateRWF(celli);
    }
    
    forAll(RWF_.boundaryField(), patchi)
    {
        fvPatchScalarField& pRWF = RWF_.boundaryFieldRef()[patchi];
        
        forAll(pRWF, facei)
        {
            pRWF[facei] = recalculatepRWF(patchi, facei);
        }
    }
}


void dsmcSpherical::writeSphericalInfo() const
{
    Info<< nl << "Spherical simulation:" << nl
        << "- coordinate system origin" << tab << origin_ << nl
        << "- radial weighting method" << tab << rWMethod_ << "-based" << nl
        << "- radial extent" << tab << radialExtent_ << nl
        << "- maximum radial weighting factor" << tab << maxRWF_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcSpherical::dsmcSpherical
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcCoordinateSystem(t, mesh, cloud),
    cloud_(cloud),
    radialExtent_(0.0),
    maxRWF_(1.0),
    origin_(vector::zero),
    rWMethod_(word::null),
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcSpherical::~dsmcSpherical()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcSpherical::checkCoordinateSystemInputs(const bool init)
{
    rWMethod_ = cloud_.particleProperties().subDict("sphericalProperties")
        .lookupOrDefault<word>("radialWeightingMethod", "cell");
    
    if
    (
        rWMethod_ != "cell" and rWMethod_ != "particle"
            and rWMethod_ != "mixed"
    )
    {
        FatalErrorIn
        (
            "dsmcSpherical::checkCoordinateSystemInputs(const bool init)"
        )
        << "The radial weighting method is badly defined. Choices "
           "in constant/dsmcProperties are cell, particle, or "
           "mixed. Please edit the entry: radialWeightingMethod"
        << exit(FatalError);
    }
    
    maxRWF_ = readScalar
        (
            cloud_.particleProperties().subDict("sphericalProperties")
                .lookup("maxRadialWeightingFactor")
        );
        
    origin_ = cloud_.particleProperties().subDict("sphericalProperties")
        .lookupOrDefault<vector>("origin", vector::zero);
        
    scalarField radii(mesh_.faceCentres().size(), 0.0);
    
    forAll(mesh_.faceCentres(), i)
    {
        radii[i] = sqrt
            (
                sqr(mesh_.faceCentres()[i].x() - origin_.x())
              + sqr(mesh_.faceCentres()[i].y() - origin_.y())
              + sqr(mesh_.faceCentres()[i].z() - origin_.z())
            );
    }
        
    radialExtent_ = gMax(radii);
        
    writeCoordinateSystemInfo();
         
    if (init)
    {
        // "particle" cannot be used in dsmcInitialise, "cell" is thus employed
        rWMethod_ = "cell"; 
    }
    else
    {
        if (rWMethod_ != "particle")
        {
            updateRWF();
        }
    }
}


void dsmcSpherical::evolve()
{
    if (rWMethod_ == "particle")
    {
        updateRWF();
    }
    
    sphericalWeighting();
    cloud_.reBuildCellOccupancy();
}


scalar dsmcSpherical::recalculatepRWF
(
    const label patchI,
    const label faceI
) const
{
    const point& fC = mesh_.boundaryMesh()[patchI].faceCentres()[faceI];
    const scalar radius = 
        sqrt
        (
            sqr(fC.x() - origin_.x())
          + sqr(fC.y() - origin_.y())
          + sqr(fC.z() - origin_.z())
        );
    
    return 1.0 + (maxRWF() - 1.0)*sqr(radius/radialExtent()); 
}


scalar dsmcSpherical::recalculateRWF
(
    const label cellI, 
    const bool mixedRWMethod
) const
{
    scalar RWF = 1.0;
    
    if (rWMethod_ == "particle" or (mixedRWMethod and rWMethod_ == "mixed"))
    {
        const DynamicList<dsmcParcel*>& cellParcels(cloud_.cellOccupancy()[cellI]);
        
        RWF = 0.0;
        label nMols = 0;
        
        forAll(cellParcels, i)
        {
            const dsmcParcel& p = *cellParcels[i];
            
            const scalar radius = 
                sqrt
                (
                    sqr(p.position().x() - origin_.x())
                  + sqr(p.position().y() - origin_.y())
                  + sqr(p.position().z() - origin_.z())
                );

            RWF += 1.0 + (maxRWF() - 1.0)*sqr(radius/radialExtent());
            
            nMols++;
        }
        
        RWF /= max(nMols, 1);
    }
    else
    {
        const point& cC = mesh_.cellCentres()[cellI];
        const scalar radius = 
            sqrt
            (
                sqr(cC.x() - origin_.x())
              + sqr(cC.y() - origin_.y())
              + sqr(cC.z() - origin_.z())
            );
    
        RWF += (maxRWF() - 1.0)*sqr(radius/radialExtent());
    }
    
    return RWF;    
}


void dsmcSpherical::writeCoordinateSystemInfo() const
{
    writeSphericalInfo();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
