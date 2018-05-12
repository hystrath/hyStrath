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
    dsmcGlobalPorousMediumMeasurements

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcGlobalPorousMediumMeasurements.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcGlobalPorousMediumMeasurements, 0);

    addToRunTimeSelectionTable
    (
        porousMeasurements, 
        dsmcGlobalPorousMediumMeasurements,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcGlobalPorousMediumMeasurements::updateMediumPropertiesMeasurement
(
    dsmcParcel& p,
    const label& delPatchId
)
{
    if (p.isTracked())
    {
        if (p.tracked().inPatchId() != -1)
        {
            p.tracked().updateTotalDistanceTravelled(p.position());
            
            const word patchesKey = std::to_string(p.tracked().inPatchId()) 
                + std::to_string(delPatchId);
            
            if (tracerPatchNamesMap_.found(patchesKey))
            {
                //- The deletion patch is in the list of outflow patches of the
                //  tracer inflow patch
                dsmcParcel::TrackedParcel::nDELETED++;
                
                mediumTotalDistanceTravelled_ += 
                    p.tracked().distanceTravelled();
                
                mediumTransitTime_ += mesh_.time().value()
                    - p.tracked().initialTime();
            }
            else if (patchesKey[0] == patchesKey[1])
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


void dsmcGlobalPorousMediumMeasurements::updateMediumPropertiesMeasurement_cyclic
(
    dsmcParcel& p,
    const label& neiPatchId,
    const label& orgPatchId,
    const vector& orgPosition
)
{
    /*if (measureMediumProperties_)
    {
        if (p.isTracked())
        {
            // orgPatchId is the id of the patch on which the particle strikes, i.e., the 'deletion' patch
            // neiPatchId is the id of the neighbouring patch, i.e., the patch of reintroduction
            
            // TODO update mediumTotalDistanceTravelled_
            
            if (orgPatchId != p.tracked().inPatchId())
            {
                //- The cyclic patch is different from the patch of insertion
                //Info << "inP: " << p.tracked().inPatchId() << tab << "orgP: " << orgPatchId << tab << "neiP: " << neiPatchId << endl;
                dsmcParcel::TrackedParcel::nDELETED++;
                
                mediumTotalDistanceTravelled_ += p.tracked().distanceTravelled();
                    
                mediumTransitTime_ += mesh_.time().value() - p.tracked().initialTime();
                
                if (trackingProbability_ > rndGen_.sample01<scalar>())
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
                    if (orgPatchId == mesh_.boundaryMesh().findPatchID(trackFromPatchNames_[patchNameCounter]))
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
            if (trackingProbability_ > rndGen_.sample01<scalar>())
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


void dsmcGlobalPorousMediumMeasurements::writeGlobalPorousMeasurementsInfo() const
{
    Info<< nl << "Porous media measurements - global:" << nl
        << "- tracking probability" << tab << trackingProbability_ << nl
        << "- dimensionality" << tab << dimensionality_ << nl
        << "- map of tracer patch names" << tab << tracerPatchNamesMap_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcGlobalPorousMediumMeasurements::dsmcGlobalPorousMediumMeasurements
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    porousMeasurements(t, mesh, cloud),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    trackingProbability_(0),
    tracerPatchNamesMap_(),
    dimensionality_(3),
    mediumTransitTime_(0.0),
    mediumTotalDistanceTravelled_(0.0),
    mediumTortuosity_(0.0),
    timeLooping_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcGlobalPorousMediumMeasurements::~dsmcGlobalPorousMediumMeasurements()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcGlobalPorousMediumMeasurements::checkPorousMeasurementsInputs()
{
    dsmcParcel::TrackedParcel::nDELETED = 0;
    nLooping_ = 0;
    
    if (cloud_.particleProperties().isDict("tracerProperties"))
    {
        wordList trackFromPatchNames;
        List<wordList> trackToPatchNames;
        
        if
        (
            cloud_.particleProperties().subDict("tracerProperties")
              .found("inflowPatchNames")
        )
        {
            cloud_.particleProperties().subDict("tracerProperties")
                .lookup("inflowPatchNames") >> trackFromPatchNames;
        }
            
        if
        (
            cloud_.particleProperties().subDict("tracerProperties")
              .found("outflowPatchNames")
        )
        {
            cloud_.particleProperties().subDict("tracerProperties")
               .lookup("outflowPatchNames") >> trackToPatchNames;
        }
            
        
        if (trackFromPatchNames.size() != trackToPatchNames.size())
        {
            FatalErrorIn
            (
                "dsmcGlobalPorousMediumMeasurements::"
                "checkPorousMeasurementsInputs()"
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
            if (mesh_.boundaryMesh().findPatchID(trackFromPatchNames[inflowPatch]) == -1)
            {
                FatalErrorIn
                (
                    "dsmcGlobalPorousMediumMeasurements::"
                    "checkPorousMeasurementsInputs()"
                )
                << "The tracers inflow patch " << trackFromPatchNames[inflowPatch]
                << " defined in constant/dsmcProperties does not exist"
                << exit(FatalError);
            }
            
            forAll(trackToPatchNames[inflowPatch], outflowPatch)
            {
                if (mesh_.boundaryMesh().findPatchID(trackToPatchNames[inflowPatch][outflowPatch]) == -1)
                {
                    FatalErrorIn
                    (
                        "dsmcGlobalPorousMediumMeasurements::"
                        "checkPorousMeasurementsInputs()"
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
                const word key = 
                    std::to_string
                    (
                        mesh_.boundaryMesh().findPatchID
                        (
                            trackFromPatchNames[inflowPatch]
                        )
                    )
                  + std::to_string
                    (
                        mesh_.boundaryMesh().findPatchID
                        (
                            trackToPatchNames[inflowPatch][outflowPatch]
                        )
                    );
                
                tracerPatchNamesMap_.insert(key);
            }
        }
        
        trackingProbability_ = cloud_.particleProperties()
            .subDict("tracerProperties")
            .lookupOrDefault<scalar>("trackingProbability", 0.1);
            
        dimensionality_ = cloud_.particleProperties()
            .subDict("tracerProperties")
            .lookupOrDefault<label>("dimensionality", 3);
    }
    else
    {
        FatalErrorIn
        (
            "dsmcGlobalPorousMediumMeasurements::"
            "checkPorousMeasurementsInputs()"
        )
        << "Subdictionary tracerProperties is missing in "
           "constant/dsmcProperties"
        << endl;
    }
        
    writeGlobalPorousMeasurementsInfo();
}


void dsmcGlobalPorousMediumMeasurements::update()
{}


void dsmcGlobalPorousMediumMeasurements::specularInteraction
(
    dsmcParcel& p,
    const vector& nw
)
{
    if (p.isTracked())
    {
        if (p.tracked().inPatchId() != -1)
        {
            p.tracked().updateTotalDistanceTravelled(p.position());
        }
    }
}


void dsmcGlobalPorousMediumMeasurements::diffuseInteraction(dsmcParcel& p)
{
    if (p.isTracked())
    {
        if (p.tracked().inPatchId() != -1)
        {
            p.tracked().updateTotalDistanceTravelled(p.position());
        }
    }
}


void dsmcGlobalPorousMediumMeasurements::deletionInteraction
(
    dsmcParcel& p,
    const label patchId
)
{
    updateMediumPropertiesMeasurement(p, patchId);
}


void dsmcGlobalPorousMediumMeasurements::additionInteraction
(
    dsmcParcel& p,
    const label patchId
)
{
    if (patchId != -1)
    {
        if (trackingProbability_ > cloud_.rndGen().sample01<scalar>())
        {
            p.setTracked
            (
                true, 
                patchId,
                mesh_.time().value(), 
                p.position()
            );
        }
    }
}


void dsmcGlobalPorousMediumMeasurements::writePorousMeasurementsInfo() const
{
    const scalar nDel = dsmcParcel::TrackedParcel::nDELETED + SMALL;
    const scalar nLooping = nLooping_ + SMALL;
    
    scalar weightingFactor = 0.0;
    
    if (dsmcParcel::TrackedParcel::nDELETED > 0)
    {
        weightingFactor = nLooping_/nDel;
    }
    
    /*const scalar mediumSpatialExtension = 1.0; // TODO
    
    mediumTortuosity_ = mediumTotalDistanceTravelled_
        /(nDel*mediumSpatialExtension);*/
    
    Info << "    Tracked parcels recorded        = " 
         << dsmcParcel::TrackedParcel::nDELETED << nl
         << "    Mean particle displacement      = " 
         << mediumTotalDistanceTravelled_/nDel << nl
         //<< "    Medium tortuosity               = " 
         //<< mediumTortuosity_ << nl
         << "    Recorded flow transit time      = "
         << mediumTransitTime_/nDel << nl
         << "    Looping parcels recorded        = " 
         << nLooping_ << nl
         << "    Flow time penalty               = "
         << timeLooping_/nLooping << nl 
         << "    Corrected flow transit time     = "
         << mediumTransitTime_/nDel + weightingFactor*timeLooping_/nLooping
         << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
