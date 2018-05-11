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
    dsmcAxisymmetric

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcAxisymmetric.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcAxisymmetric, 0);

    addToRunTimeSelectionTable
    (
        dsmcCoordinateSystem, 
        dsmcAxisymmetric,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcAxisymmetric::axisymmetricWeighting()
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
                    
                    U.component(angularCoordinate_) *= -1.0;
                    
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
                    
                    U.component(angularCoordinate_) *= -1.0;

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


void dsmcAxisymmetric::updateRWF()
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


void dsmcAxisymmetric::writeAxisymmetricInfo() const
{
    Info<< nl << "Axisymmetric simulation:" << nl
        << "- revolution axis label" << tab << revolutionAxis_ << nl
        << "- polar axis label" << tab << polarAxis_ << nl
        << "- angular coordinate label" << tab << angularCoordinate_ << nl 
        << "- radial weighting method" << tab << rWMethod_ << "-based" << nl
        << "- radial extent" << tab << radialExtent_ << nl
        << "- maximum radial weighting factor" << tab << maxRWF_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcAxisymmetric::dsmcAxisymmetric
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcCoordinateSystem(t, mesh, cloud),
    cloud_(cloud),
    revolutionAxis_(0),
    polarAxis_(1),
    angularCoordinate_(2),
    radialExtent_(0.0),
    maxRWF_(1.0),
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

dsmcAxisymmetric::~dsmcAxisymmetric()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAxisymmetric::checkCoordinateSystemInputs(const bool init)
{
    rWMethod_ = cloud_.particleProperties().subDict("axisymmetricProperties")
        .lookupOrDefault<word>("radialWeightingMethod", "cell");
    
    if
    (
        rWMethod_ != "cell" and rWMethod_ != "particle"
            and rWMethod_ != "mixed"
    )
    {
        FatalErrorIn
        (
            "dsmcAxisymmetric::checkCoordinateSystemInputs(const bool init)"
        )
        << "The radial weighting method is badly defined. Choices "
           "in constant/dsmcProperties are cell, particle, or "
           "mixed. Please edit the entry: radialWeightingMethod"
        << exit(FatalError);
    }
    
    const word& revolutionAxis = 
        cloud_.particleProperties().subDict("axisymmetricProperties")
            .lookupOrDefault<word>("revolutionAxis", word::null);
        
    const word& polarAxis = 
        cloud_.particleProperties().subDict("axisymmetricProperties")
            .lookupOrDefault<word>("polarAxis", word::null);
        
    if (revolutionAxis == "z")
    {
        revolutionAxis_ = 2;
        
        if (polarAxis == word::null) 
        {
            polarAxis_ = (revolutionAxis_ + 1)%3;
            angularCoordinate_ = (revolutionAxis_ + 2)%3;
        }
        else if (polarAxis == "y")
        {
            polarAxis_ = 1;
            angularCoordinate_ = 0;
        }
        else if (polarAxis == "x")
        {
            polarAxis_ = 0;
            angularCoordinate_ = 1;
        }
        else
        {
            FatalErrorIn
            (
                "dsmcAxisymmetric::checkCoordinateSystemInputs(const bool init)"
            )
            << "Revolution and polar axes are badly defined in "
               "constant/dsmcProperties axisymmetricProperties{}"
            << exit(FatalError);
        }
    }
    else if (revolutionAxis == "y")
    {
        revolutionAxis_ = 1;
        
        if (polarAxis == word::null) 
        {
            polarAxis_ = (revolutionAxis_ + 1)%3;
            angularCoordinate_ = (revolutionAxis_ + 2)%3;
        }
        else if (polarAxis == "x")
        {
            polarAxis_ = 0;
            angularCoordinate_ = 2;
        }
        else if (polarAxis == "z")
        {
            polarAxis_ = 2;
            angularCoordinate_ = 0;
        }
        else
        {
            FatalErrorIn
            (
                "dsmcAxisymmetric::checkCoordinateSystemInputs(const bool init)"
            )
            << "Revolution and polar axes are badly defined in "
               "constant/dsmcProperties axisymmetricProperties{}"
            << exit(FatalError);
        }
    }
    else if (revolutionAxis == "x")
    {
        if (polarAxis == "z")
        {
            polarAxis_ = 2;
            angularCoordinate_ = 1;
        }
        else if (polarAxis != "y")
        {
            FatalErrorIn
            (
                "dsmcAxisymmetric::checkCoordinateSystemInputs(const bool init)"
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
        
    if (!(radialExtent_ > 0))
    {
        radialExtent_ = -gMin
        (
            mesh_.faceCentres().component(polarAxis_)
        );
    }
        
    maxRWF_ = readScalar
        (
            cloud_.particleProperties().subDict("axisymmetricProperties")
                .lookup("maxRadialWeightingFactor")
        );
    
    writeCoordinateSystemInfo();
         
    if (init)
    {
        // "particle" cannot be used in dsmcInitialise, "cell" is thus employed
        rWMethod_ = "cell"; 
        
        updateRWF();
    }
    else
    {
        if (rWMethod_ != "particle")
        {
            updateRWF();
        }
    }
}


void dsmcAxisymmetric::evolve()
{
    if (rWMethod_ == "particle")
    {
        updateRWF();
    }
    
    axisymmetricWeighting();
    cloud_.reBuildCellOccupancy();
}


scalar dsmcAxisymmetric::recalculatepRWF
(
    const label patchI,
    const label faceI
) const
{
    const point& fC = mesh_.boundaryMesh()[patchI].faceCentres()[faceI];
    const scalar radius = mag(fC.component(polarAxis_));
    
    return 1.0 + (maxRWF() - 1.0)*radius/radialExtent(); 
}


scalar dsmcAxisymmetric::recalculateRWF
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
    
    return RWF;    
}


void dsmcAxisymmetric::writeCoordinateSystemInfo() const
{
    writeAxisymmetricInfo();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
