/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extractEulerianParticles.H"
#include "regionSplit2D.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "pairPatchAgglomeration.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "binned.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(extractEulerianParticles, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        extractEulerianParticles,
        dictionary
    );
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileName
Foam::functionObjects::extractEulerianParticles::dictBaseFileDir() const
{
    fileName baseDir(".."); //  = mesh_.time().path();

    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        baseDir = baseDir/".."/"postProcessing";
    }
    else
    {
        baseDir = baseDir/"postProcessing";
    }

    return baseDir;
}


void Foam::functionObjects::extractEulerianParticles::checkFaceZone()
{
    DebugInFunction << endl;

    zoneID_ = mesh_.faceZones().findZoneID(faceZoneName_);
    if (zoneID_ == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName_
            << ".  Available faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const label nFaces = fz.size();
    const label allFaces = returnReduce(nFaces, sumOp<label>());

    if (allFaces < nInjectorLocations_)
    {
        FatalErrorInFunction
            << "faceZone " << faceZoneName_
            << ": Number of faceZone faces (" << allFaces
            << ") is less than the number of requested locations ("
            << nInjectorLocations_ << ")."
            << exit(FatalError);
    }

    Info<< type() << " " << name() << " output:" << nl
        << "    faceZone : " << faceZoneName_ << nl
        << "    faces    : " << allFaces << nl
        << endl;

    // Initialise old iteration blocked faces
    // Note: for restart, this info needs to be written/read
    regions0_.setSize(fz.size(), -1);
}


void Foam::functionObjects::extractEulerianParticles::initialiseBins()
{
    DebugInFunction << endl;

    if (!nInjectorLocations_)
    {
        return;
    }

    const faceZone& fz = mesh_.faceZones()[zoneID_];

    // Agglomerate faceZone faces into nInjectorLocations_ global locations
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    const label nFaces = fz.size();
    label nLocations = nInjectorLocations_;

    if (Pstream::parRun())
    {
        label nGlobalFaces = returnReduce(nFaces, sumOp<label>());
        scalar fraction = scalar(nFaces)/scalar(nGlobalFaces);
        nLocations = ceil(fraction*nInjectorLocations_);
        if (debug)
        {
            Pout<< "nFaces:" << nFaces
                << ", nGlobalFaces:" << nGlobalFaces
                << ", fraction:" << fraction
                << ", nLocations:" << nLocations
                << endl;
        }
    }

    pairPatchAgglomeration ppa
    (
        patch.localFaces(),
        patch.localPoints(),
        10,
        50,
        nLocations,
        labelMax,
        180
    );

    ppa.agglomerate();

    label nCoarseFaces = 0;
    if (nFaces != 0)
    {
        fineToCoarseAddr_ = ppa.restrictTopBottomAddressing();
        nCoarseFaces = max(fineToCoarseAddr_) + 1;
    }

    globalCoarseFaces_ = globalIndex(nCoarseFaces);

    Info<< "Created " << returnReduce(nCoarseFaces, sumOp<label>())
        << " coarse faces" << endl;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::extractEulerianParticles::phiU() const
{
    DebugInFunction << endl;

    const surfaceScalarField& phi
    (
        mesh_.lookupObject<surfaceScalarField>(phiName_)
    );

    if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        return phi/fvc::interpolate(rho);
    }

    return phi;
}


void Foam::functionObjects::extractEulerianParticles::setBlockedFaces
(
    const surfaceScalarField& alphaf,
    const faceZone& fz,
    boolList& blockedFaces
)
{
    DebugInFunction << endl;

    // Initialise storage for patch and patch-face indices where faceZone
    // intersects mesh patch(es)
    patchIDs_.setSize(fz.size(), -1);
    patchFaceIDs_.setSize(fz.size(), -1);

    label nBlockedFaces = 0;
    forAll(fz, localFacei)
    {
        const label facei = fz[localFacei];

        if (mesh_.isInternalFace(facei))
        {
            if (alphaf[facei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }
        }
        else
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
            label patchFacei = -1;

            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            const scalarField& alphafp = alphaf.boundaryField()[patchi];

            if (isA<coupledPolyPatch>(pp))
            {
                const coupledPolyPatch& cpp =
                    refCast<const coupledPolyPatch>(pp);

                if (cpp.owner())
                {
                    patchFacei = cpp.whichFace(facei);
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                patchFacei = pp.whichFace(facei);
            }

            if (patchFacei == -1)
            {
                patchi = -1;
            }
            else if (alphafp[patchFacei] > alphaThreshold_)
            {
                blockedFaces[localFacei] = true;
            }

            patchIDs_[localFacei] = patchi;
            patchFaceIDs_[localFacei] = patchFacei;
        }
    }

    DebugInFunction << "Number of blocked faces: " << nBlockedFaces << endl;
}


void Foam::functionObjects::extractEulerianParticles::collectParticles
(
    const scalar time,
    const boolList& collectParticle
)
{
    DebugInFunction << "collectParticle: " << collectParticle << endl;

    // Collect particles on local processors that have passed through faceZone
    forAll(collectParticle, regioni)
    {
        if (!collectParticle[regioni])
        {
            continue;
        }

        Map<label>::const_iterator iter = regionToParticleMap_.find(regioni);
        eulerianParticle p = particles_[iter()];

        if (p.faceIHit != -1 && nInjectorLocations_)
        {
            // Use coarse face index for tag output
            label coarseFacei = fineToCoarseAddr_[p.faceIHit];
            p.faceIHit = globalCoarseFaces_.toGlobal(coarseFacei);
        }

        reduce(p, sumParticleOp<eulerianParticle>());

        const scalar pDiameter = cbrt(6.0*p.V/constant::mathematical::pi);

        if ((pDiameter > minDiameter_) && (pDiameter < maxDiameter_))
        {
            if (Pstream::master())
            {
                const scalar d = cbrt(6*p.V/constant::mathematical::pi);
                const point position = p.VC/(p.V + ROOTVSMALL);
                const vector U = p.VU/(p.V + ROOTVSMALL);
                label tag = -1;
                if (nInjectorLocations_)
                {
                    tag = p.faceIHit;
                }

                injectedParticle* ip = new injectedParticle
                (
                    mesh_,
                    position,
                    tag,
                    time,
                    d,
                    U
                );

                cloud_.addParticle(ip);
            }

            nCollectedParticles_++;
        }
        else
        {
            // Discard particles over/under diameter threshold
            nDiscardedParticles_++;
            discardedVolume_ += p.V;
        }
    }
}


void Foam::functionObjects::extractEulerianParticles::calculateAddressing
(
    const label nRegionsOld,
    const label nRegionsNew,
    const scalar time,
    labelList& regionFaceIDs
)
{
    DebugInFunction << endl;

    // New region can only point to one old region
    // Old region can only point to one new region.  If old region intersects
    // multiple new regions, select max of new region indices.

    labelList oldToNewRegion(nRegionsOld, -1);
    labelList newToOldRegion(nRegionsNew, -1);

    forAll(regionFaceIDs, facei)
    {
        label newRegioni = regionFaceIDs[facei];
        label oldRegioni = regions0_[facei];

        if (newRegioni != -1)
        {
            newToOldRegion[newRegioni] = oldRegioni;

            if (oldRegioni != -1)
            {
                // New region linked to old (so can accumulate particle data)
                // Old region might already be mapped to a new region
                oldToNewRegion[oldRegioni] =
                    max(newRegioni, oldToNewRegion[oldRegioni]);
            }
        }
    }

    // Need to re-number the new region indices based on knowledge of which
    // old region they intersect.  After operations, there should be a
    // one-to-one match between the old and new regions.

    // Ensure all old regions point to the same new regions (if present)
    Pstream::listCombineGather(oldToNewRegion, maxEqOp<label>());
    Pstream::listCombineScatter(oldToNewRegion);

    // Any new region that points to an old region should be renumbered to the
    // new region specified by the oldToNewRegion index

    if (oldToNewRegion.size())
    {
        // Create corrected new to new addressing
        labelList newToNewRegionCorr(newToOldRegion.size(), -1);
        forAll(newToOldRegion, newRegioni)
        {
            label oldRegioni = newToOldRegion[newRegioni];
            if (oldRegioni != -1)
            {
                label newRegionICorr = oldToNewRegion[oldRegioni];
                newToNewRegionCorr[newRegioni] = newRegionICorr;
                newToOldRegion[newRegioni] = oldRegioni;
            }
        }

        // Renumber the new (current) face region IDs
        forAll(regionFaceIDs, facei)
        {
            label newRegioni = regionFaceIDs[facei];

            if (newRegioni != -1)
            {
                label newRegionICorr = newToNewRegionCorr[newRegioni];

                // Note: newRegionICorr can be -1 if the region is new, since
                // the address corrections are based on inverting the
                // old-to-new addressing
                if (newRegionICorr != -1)
                {
                    regionFaceIDs[facei] = newRegionICorr;
                }
            }
        }


        boolList collectParticleFlag(nRegionsOld, true);
        forAll(oldToNewRegion, oldRegioni)
        {
            label newRegioni = oldToNewRegion[oldRegioni];
            if (newRegioni != -1)
            {
                collectParticleFlag[oldRegioni] = false;
            }
        }

        // Collect particles whose IDs are no longer active
        collectParticles(time, collectParticleFlag);
    }


    // Re-order collection bins

    Map<label> newRegionToParticleMap(nRegionsNew);
    List<eulerianParticle> newParticles(nRegionsNew);
    label particlei = 0;

    forAll(newToOldRegion, newRegioni)
    {
        label oldRegioni = newToOldRegion[newRegioni];
        if (oldRegioni == -1)
        {
            // No mapping from old to new - this is a new particle
            newRegionToParticleMap.insert(newRegioni, particlei);
            particlei++;
        }
        else
        {
            // Update addressing for old to new regions
            label oldParticlei = regionToParticleMap_[oldRegioni];
            if (newRegionToParticleMap.insert(newRegioni, particlei))
            {
                // First time this particle has been seen
                newParticles[particlei] = particles_[oldParticlei];
                particlei++;
            }
            else
            {
                // Combine with existing contribution
                label newParticlei = newRegionToParticleMap[newRegioni];
                newParticles[newParticlei] =
                    sumParticleOp<eulerianParticle>()
                    (
                        newParticles[newParticlei],
                        particles_[oldParticlei]
                    );
            }
        }
    }

    // Reset the particles list and addressing for latest available info
    particles_.transfer(newParticles);
    regionToParticleMap_ = newRegionToParticleMap;
}


void Foam::functionObjects::extractEulerianParticles::accumulateParticleInfo
(
    const surfaceScalarField& alphaf,
    const surfaceScalarField& phi,
    const labelList& regionFaceIDs,
    const faceZone& fz
)
{
    DebugInFunction << endl;

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const surfaceVectorField Uf(fvc::interpolate(U));

    const scalar deltaT = mesh_.time().deltaTValue();
    const pointField& faceCentres = mesh_.faceCentres();

    forAll(regionFaceIDs, localFacei)
    {
        const label newRegioni = regionFaceIDs[localFacei];

        if (newRegioni != -1)
        {
            const label particlei = regionToParticleMap_[newRegioni];
            const label meshFacei = fz[localFacei];
            eulerianParticle& p = particles_[particlei];

            if (p.faceIHit < 0)
            {
                // New particle - does not exist in particles_ list
                p.faceIHit = localFacei;
                p.V = 0;
                p.VC = vector::zero;
                p.VU = vector::zero;
            }

            // Accumulate particle properties
            scalar magPhii = mag(faceValue(phi, localFacei, meshFacei));
            vector Ufi = faceValue(Uf, localFacei, meshFacei);
            scalar dV = magPhii*deltaT;
            p.V += dV;
            p.VC += dV*faceCentres[meshFacei];
            p.VU += dV*Ufi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticles::extractEulerianParticles
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name),
    cloud_(mesh_, "eulerianParticleCloud"),
    faceZoneName_(word::null),
    zoneID_(-1),
    patchIDs_(),
    patchFaceIDs_(),
    alphaName_("alpha"),
    alphaThreshold_(0.1),
    UName_("U"),
    rhoName_("rho"),
    phiName_("phi"),
    nInjectorLocations_(0),
    fineToCoarseAddr_(),
    globalCoarseFaces_(),
    regions0_(),
    nRegions0_(0),
    particles_(),
    regionToParticleMap_(),
    minDiameter_(ROOTVSMALL),
    maxDiameter_(GREAT),
    nCollectedParticles_(0),
    nDiscardedParticles_(0),
    discardedVolume_(0)
{
    if (mesh_.nSolutionD() != 3)
    {
        FatalErrorInFunction
            << name << " function object only applicable to 3-D cases"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::extractEulerianParticles::~extractEulerianParticles()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::extractEulerianParticles::read
(
    const dictionary& dict
)
{
    DebugInFunction << endl;

    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.lookup("faceZone") >> faceZoneName_;
        dict.lookup("alpha") >> alphaName_;

        dict.readIfPresent("alphaThreshold", alphaThreshold_);
        dict.readIfPresent("U", UName_);
        dict.readIfPresent("rho", rhoName_);
        dict.readIfPresent("phi", phiName_);
        dict.readIfPresent("nLocations", nInjectorLocations_);
        dict.readIfPresent("minDiameter", minDiameter_);
        dict.readIfPresent("maxDiameter", maxDiameter_);

        checkFaceZone();

        if (nInjectorLocations_)
        {
            initialiseBins();
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::extractEulerianParticles::execute()
{
    DebugInFunction << endl;

    Log << type() << " " << name() << " output:" << nl;

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    const surfaceScalarField alphaf
    (
        typeName + ":alphaf",
        fvc::interpolate(alpha)
    );

    const faceZone& fz = mesh_.faceZones()[zoneID_];
    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), fz),
        mesh_.points()
    );

    // Set the blocked faces, i.e. where alpha > alpha threshold value
    boolList blockedFaces(fz.size(), false);
    setBlockedFaces(alphaf, fz, blockedFaces);

    // Split the  faceZone according to the blockedFaces
    // - Returns a list of (disconnected) region index per face zone face
    regionSplit2D regionFaceIDs(mesh_, patch, blockedFaces);

    // Global number of regions
    const label nRegionsNew = regionFaceIDs.nRegions();

    // Calculate the addressing between the old and new region information
    // Also collects particles that have traversed the faceZone
    calculateAddressing
    (
        nRegions0_,
        nRegionsNew,
        mesh_.time().value(),
        regionFaceIDs
    );

    // Process latest region information
    tmp<surfaceScalarField> tphi = phiU();
    accumulateParticleInfo(alphaf, tphi(), regionFaceIDs, fz);

    // Reset the blocked faces for the next integration step
    nRegions0_ = nRegionsNew;
    regions0_ = regionFaceIDs;

    Log << "    Collected particles: " << nCollectedParticles_ << nl
        << "    Discarded particles: " << nDiscardedParticles_ << nl
        << "    Discarded volume: " << discardedVolume_ << nl
        << endl;

    return true;
}


bool Foam::functionObjects::extractEulerianParticles::write()
{
    DebugInFunction << endl;

    cloud_.write();

    return true;
}


// ************************************************************************* //
