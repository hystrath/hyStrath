/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "fluxSummary.H"
#include "surfaceFields.H"
#include "surfFields.H"
#include "surfMesh.H"
#include "dictionary.H"
#include "Time.H"
#include "syncTools.H"
#include "meshTools.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"
#include "globalIndex.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluxSummary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        fluxSummary,
        dictionary
    );
}
template<>
const char* NamedEnum
<
    functionObjects::fluxSummary::modeType,
    5
>::names[] =
{
    "faceZone",
    "faceZoneAndDirection",
    "cellZoneAndDirection",
    "surface",
    "surfaceAndDirection"
};
}


const Foam::NamedEnum<Foam::functionObjects::fluxSummary::modeType, 5>
Foam::functionObjects::fluxSummary::modeTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fluxSummary::isSurfaceMode() const
{
    bool isSurf = false;

    switch (mode_)
    {
        case mdSurface:
        case mdSurfaceAndDirection:
            isSurf = true;
            break;

        default:
            break;
    }

    return isSurf;
}


Foam::word Foam::functionObjects::fluxSummary::checkFlowType
(
    const dimensionSet& fieldDims,
    const word& fieldName
) const
{
    // Surfaces are multipled by their area, so account for that
    // in the dimension checking
    dimensionSet dims =
        fieldDims * (isSurfaceMode() ? dimTime*dimArea : dimTime);

    if (dims == dimVolume)
    {
        return "volumetric";
    }
    else if (dims == dimMass)
    {
        return "mass";
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported flux field " << fieldName << " with dimensions "
            << fieldDims
            << ".  Expected either mass flow or volumetric flow rate."
            << abort(FatalError);

        return Foam::word::null;
    }
}


void Foam::functionObjects::fluxSummary::initialiseSurface
(
    const word& surfName,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<boolList>& faceFlip
) const
{
    const surfMesh* sPtr = mesh_.lookupObjectPtr<surfMesh>(surfName);
    if (!sPtr)
    {
        FatalErrorInFunction
            << "Unable to find surface " << surfName
            << ".  Valid surfaces are: " << mesh_.sortedNames<surfMesh>()
            << '.'
            << exit(FatalError);
    }

    names.append(surfName);
    directions.append(Zero); // dummy value
    faceFlip.append(boolList(0)); // no flip-map
}


void Foam::functionObjects::fluxSummary::initialiseSurfaceAndDirection
(
    const word& surfName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<boolList>& faceFlip
) const
{
    const surfMesh* sPtr = mesh_.lookupObjectPtr<surfMesh>(surfName);
    if (!sPtr)
    {
        FatalErrorInFunction
            << "Unable to find surface " << surfName
            << ".  Valid surfaces are: " << mesh_.sortedNames<surfMesh>()
            << '.'
            << exit(FatalError);
    }

    const surfMesh& s = *sPtr;
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    names.append(surfName);
    directions.append(refDir);
    faceFlip.append(boolList(0));

    boolList& flips = faceFlip[faceFlip.size()-1];
    flips.setSize(s.size(), false);

    forAll(s, i)
    {
        // orientation set by comparison with reference direction
        const vector& n = s.faceNormals()[i];

        if ((n & refDir) > tolerance_)
        {
            flips[i] = false;
        }
        else
        {
            flips[i] = true;
        }
    }
}


void Foam::functionObjects::fluxSummary::initialiseFaceZone
(
    const word& faceZoneName,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    label zonei = mesh_.faceZones().findZoneID(faceZoneName);
    if (zonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }
    const faceZone& fZone = mesh_.faceZones()[zonei];

    names.append(faceZoneName);
    directions.append(Zero); // dummy value

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<bool>  flips(fZone.size());

    forAll(fZone, i)
    {
        label facei = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceID = facei;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(facei);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = facei - pp.start();
            }
            else
            {
                faceID = -1;
                facePatchID = -1;
            }
        }

        if (faceID >= 0)
        {
            // Orientation set by faceZone flip map
            if (fZone.flipMap()[facei])
            {
                flips.append(true);
            }
            else
            {
                flips.append(false);
            }

            faceIDs.append(faceID);
            facePatchIDs.append(facePatchID);
        }
    }

    // could reduce some copying here
    faceID.append(faceIDs);
    facePatchID.append(facePatchIDs);
    faceFlip.append(flips);
}


void Foam::functionObjects::fluxSummary::initialiseFaceZoneAndDirection
(
    const word& faceZoneName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    label zonei = mesh_.faceZones().findZoneID(faceZoneName);
    if (zonei == -1)
    {
         FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid faceZones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }
    const faceZone& fZone = mesh_.faceZones()[zonei];

    names.append(faceZoneName);
    directions.append(refDir);

    DynamicList<label> faceIDs(fZone.size());
    DynamicList<label> facePatchIDs(fZone.size());
    DynamicList<bool>  flips(fZone.size());

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    vector n(Zero);

    forAll(fZone, i)
    {
        label facei = fZone[i];

        label faceID = -1;
        label facePatchID = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceID = facei;
            facePatchID = -1;
        }
        else
        {
            facePatchID = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchID];
            if (isA<coupledPolyPatch>(pp))
            {
                if (refCast<const coupledPolyPatch>(pp).owner())
                {
                    faceID = pp.whichFace(facei);
                }
                else
                {
                    faceID = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceID = facei - pp.start();
            }
            else
            {
                faceID = -1;
                facePatchID = -1;
            }
        }

        if (faceID >= 0)
        {
            // orientation set by comparison with reference direction
            if (facePatchID != -1)
            {
                n = Sf.boundaryField()[facePatchID][faceID]
                   /(magSf.boundaryField()[facePatchID][faceID] + ROOTVSMALL);
            }
            else
            {
                n = Sf[faceID]/(magSf[faceID] + ROOTVSMALL);
            }

            if ((n & refDir) > tolerance_)
            {
                flips.append(false);
            }
            else
            {
                flips.append(true);
            }

            faceIDs.append(faceID);
            facePatchIDs.append(facePatchID);
        }
    }

    // could reduce copying here
    faceID.append(faceIDs);
    facePatchID.append(facePatchIDs);
    faceFlip.append(flips);
}


void Foam::functionObjects::fluxSummary::initialiseCellZoneAndDirection
(
    const word& cellZoneName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<List<label>>& faceID,
    DynamicList<List<label>>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    const label cellZonei = mesh_.cellZones().findZoneID(cellZoneName);
    if (cellZonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find cellZone " << cellZoneName
            << ". Valid zones are: " << mesh_.cellZones().names()
            << exit(FatalError);
    }

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    labelList cellAddr(mesh_.nCells(), -1);
    const labelList& cellIDs = mesh_.cellZones()[cellZonei];
    UIndirectList<label>(cellAddr, cellIDs) = identity(cellIDs.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct boundary values for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Collect faces
    DynamicList<label> faceIDs(floor(0.1*mesh_.nFaces()));
    DynamicList<label> facePatchIDs(faceIDs.size());
    DynamicList<label> faceLocalPatchIDs(faceIDs.size());
    DynamicList<bool>  flips(faceIDs.size());

    // Internal faces
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if (((own != -1) && (nbr == -1)) || ((own == -1) && (nbr != -1)))
        {
            vector n = mesh_.faces()[facei].normal(mesh_.points());
            n /= mag(n) + ROOTVSMALL;

            if ((n & refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                flips.append(false);
            }
            else if ((n & -refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                flips.append(true);
            }
        }
    }

    // Loop over boundary faces
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, localFacei)
        {
            const label facei = pp.start() + localFacei;
            const label own = cellAddr[mesh_.faceOwner()[facei]];
            const label nbr = nbrFaceCellAddr[facei - nInternalFaces];

            if ((own != -1) && (nbr == -1))
            {
                vector n = mesh_.faces()[facei].normal(mesh_.points());
                n /= mag(n) + ROOTVSMALL;

                if ((n & refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    flips.append(false);
                }
                else if ((n & -refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    flips.append(true);
                }
            }
        }
    }

    // Convert into primitivePatch for convenience
    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), faceIDs),
        mesh_.points()
    );

    if (debug)
    {
        OBJstream os(mesh_.time().path()/"patch.obj");
        faceList faces(patch);
        os.write(faces, mesh_.points(), false);
    }


    // Data on all edges and faces
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());

    bool search = true;

    DebugInfo
        << "initialiseCellZoneAndDirection: "
        << "Starting walk to split patch into faceZones"
        << endl;

    globalIndex globalFaces(patch.size());

    label oldFaceID = 0;
    label regioni = 0;
    while (search)
    {
        DynamicList<label> changedEdges;
        DynamicList<patchEdgeFaceRegion> changedInfo;

        label seedFacei = labelMax;
        for (; oldFaceID < patch.size(); oldFaceID++)
        {
            if (allFaceInfo[oldFaceID].region() == -1)
            {
                seedFacei = globalFaces.toGlobal(oldFaceID);
                break;
            }
        }
        reduce(seedFacei, minOp<label>());

        if (seedFacei == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(seedFacei))
        {
            label localFacei = globalFaces.toLocal(seedFacei);
            const labelList& fEdges = patch.faceEdges()[localFacei];

            forAll(fEdges, i)
            {
                if (allEdgeInfo[fEdges[i]].region() != -1)
                {
                    WarningInFunction
                        << "Problem in edge face wave: attempted to assign a "
                        << "value to an edge that has already been visited. "
                        << "Edge info: " << allEdgeInfo[fEdges[i]]
                        << endl;
                }

                changedEdges.append(fEdges[i]);
                changedInfo.append(regioni);
            }
        }


        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchEdgeFaceRegion
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        if (debug)
        {
            label nCells = 0;
            forAll(allFaceInfo, facei)
            {
                if (allFaceInfo[facei].region() == regioni)
                {
                    nCells++;
                }
            }

            Info<< "*** region:" << regioni
                << "  found:" << returnReduce(nCells, sumOp<label>())
                << " faces" << endl;
        }

        regioni++;
    }

    // Collect the data per region
    label nRegion = regioni;

    List<DynamicList<label>> regionFaceIDs(nRegion);
    List<DynamicList<label>> regionFacePatchIDs(nRegion);
    List<DynamicList<bool>>  regionFaceFlips(nRegion);

    forAll(allFaceInfo, facei)
    {
        regioni = allFaceInfo[facei].region();

        regionFaceIDs[regioni].append(faceLocalPatchIDs[facei]);
        regionFacePatchIDs[regioni].append(facePatchIDs[facei]);
        regionFaceFlips[regioni].append(flips[facei]);
    }

    // Transfer to persistent storage
    forAll(regionFaceIDs, regioni)
    {
        const word zoneName = cellZoneName + ":faceZone" + Foam::name(regioni);
        names.append(zoneName);
        directions.append(refDir);
        faceID.append(regionFaceIDs[regioni]);
        facePatchID.append(regionFacePatchIDs[regioni]);
        faceFlip.append(regionFaceFlips[regioni]);

        // Write OBj of faces to file
        if (debug)
        {
            OBJstream os(mesh_.time().path()/zoneName + ".obj");
            faceList faces(mesh_.faces(), regionFaceIDs[regioni]);
            os.write(faces, mesh_.points(), false);
        }
    }

    if (log)
    {
        Info<< type() << " " << name() << " output:" << nl
            << "    Created " << faceID.size()
            << " separate face zones from cell zone " << cellZoneName << nl;

        forAll(names, i)
        {
            label nFaces = returnReduce(faceID[i].size(), sumOp<label>());
            Info<< "    " << names[i] << ": "
                << nFaces << " faces" << nl;
        }

        Info<< endl;
    }
}


Foam::scalar Foam::functionObjects::fluxSummary::totalArea
(
    const label idx
) const
{
    scalar sumMagSf = 0;

    if (isSurfaceMode())
    {
        const surfMesh& s = mesh_.lookupObject<surfMesh>(zoneNames_[idx]);
        sumMagSf = sum(s.magSf());
    }
    else
    {
        const surfaceScalarField& magSf = mesh_.magSf();

        const labelList& faceIDs = faceID_[idx];
        const labelList& facePatchIDs = facePatchID_[idx];

        forAll(faceIDs, i)
        {
            label facei = faceIDs[i];

            if (facePatchIDs[i] == -1)
            {
                sumMagSf += magSf[facei];
            }
            else
            {
                label patchi = facePatchIDs[i];
                sumMagSf += magSf.boundaryField()[patchi][facei];
            }
        }
    }

    return returnReduce(sumMagSf, sumOp<scalar>());
}


bool Foam::functionObjects::fluxSummary::surfaceModeWrite()
{
    if (zoneNames_.size())
    {
        const label surfi = 0;
        const surfMesh& s = mesh_.lookupObject<surfMesh>(zoneNames_[surfi]);
        const surfVectorField& phi = s.lookupObject<surfVectorField>(phiName_);

        Log << type() << ' ' << name() << ' '
            << checkFlowType(phi.dimensions(), phi.name()) << " write:" << nl;
    }


    forAll(zoneNames_, surfi)
    {
        const surfMesh& s = mesh_.lookupObject<surfMesh>(zoneNames_[surfi]);
        const surfVectorField& phi = s.lookupObject<surfVectorField>(phiName_);

        checkFlowType(phi.dimensions(), phi.name());

        const boolList& flips = faceFlip_[surfi];

        scalar phiPos = scalar(0);
        scalar phiNeg = scalar(0);

        tmp<scalarField> tphis = phi & s.Sf();
        const scalarField& phis = tphis();

        forAll(s, i)
        {
            scalar phif = phis[i];
            if (flips[i])
            {
                phif *= -1;
            }

            if (phif > 0)
            {
                phiPos += phif;
            }
            else
            {
                phiNeg += phif;
            }
        }

        reduce(phiPos, sumOp<scalar>());
        reduce(phiNeg, sumOp<scalar>());

        phiPos *= scaleFactor_;
        phiNeg *= scaleFactor_;

        scalar netFlux = phiPos + phiNeg;
        scalar absoluteFlux = phiPos - phiNeg;

        Log << "    surface " << zoneNames_[surfi] << ':' << nl
            << "        positive : " << phiPos << nl
            << "        negative : " << phiNeg << nl
            << "        net      : " << netFlux << nl
            << "        absolute : " << absoluteFlux
            << nl << endl;

        if (writeToFile())
        {
            filePtrs_[surfi]
                << time_.value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    Log << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fluxSummary::fluxSummary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    mode_(mdFaceZone),
    scaleFactor_(1),
    phiName_("phi"),
    zoneNames_(),
    faceID_(),
    facePatchID_(),
    faceFlip_(),
    filePtrs_(),
    tolerance_(0.8)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fluxSummary::~fluxSummary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fluxSummary::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    mode_ = modeTypeNames_.read(dict.lookup("mode"));
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    scaleFactor_ = dict.lookupOrDefault<scalar>("scaleFactor", 1.0);
    tolerance_   = dict.lookupOrDefault<scalar>("tolerance", 0.8);

    // Initialise with capacity of 10 faceZones
    DynamicList<word> faceZoneName(10);
    DynamicList<vector>       refDir(faceZoneName.capacity());
    DynamicList<List<label>>  faceID(faceZoneName.capacity());
    DynamicList<List<label>>  facePatchID(faceZoneName.capacity());
    DynamicList<boolList>     faceFlips(faceZoneName.capacity());

    switch (mode_)
    {
        case mdFaceZone:
        {
            List<word> zones(dict.lookup("faceZones"));

            forAll(zones, i)
            {
                initialiseFaceZone
                (
                    zones[i],
                    faceZoneName,
                    refDir, // fill with dummy value
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdFaceZoneAndDirection:
        {
            List<Tuple2<word, vector>>
                zoneAndDirection(dict.lookup("faceZoneAndDirection"));

            forAll(zoneAndDirection, i)
            {
                initialiseFaceZoneAndDirection
                (
                    zoneAndDirection[i].first(),
                    zoneAndDirection[i].second(),
                    faceZoneName,
                    refDir,
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdCellZoneAndDirection:
        {
            List<Tuple2<word, vector>>
                zoneAndDirection(dict.lookup("cellZoneAndDirection"));

            forAll(zoneAndDirection, i)
            {
                initialiseCellZoneAndDirection
                (
                    zoneAndDirection[i].first(),
                    zoneAndDirection[i].second(),
                    faceZoneName,
                    refDir,
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdSurface:
        {
            List<word> surfs(dict.lookup("surfaces"));

            forAll(surfs, i)
            {
                initialiseSurface
                (
                    surfs[i],
                    faceZoneName,
                    refDir,
                    faceFlips
                );
            }
            break;
        }
        case mdSurfaceAndDirection:
        {
            List<Tuple2<word, vector>>
                surfAndDirection(dict.lookup("surfaceAndDirection"));

            forAll(surfAndDirection, i)
            {
                initialiseSurfaceAndDirection
                (
                    surfAndDirection[i].first(),
                    surfAndDirection[i].second(),
                    faceZoneName,
                    refDir,
                    faceFlips
                );
            }
            break;
        }
        default:
        {
            FatalIOErrorInFunction(dict)
                << "unhandled enumeration " << modeTypeNames_[mode_]
                << abort(FatalIOError);
        }
    }

    zoneNames_.transfer(faceZoneName);
    faceID_.transfer(faceID);
    facePatchID_.transfer(facePatchID);
    faceFlip_.transfer(faceFlips);

    Info<< type() << " " << name() << " output:" << nl;

    // Calculate and report areas
    List<scalar> areas(zoneNames_.size());
    forAll(zoneNames_, zonei)
    {
        const word& zoneName = zoneNames_[zonei];
        areas[zonei] = totalArea(zonei);

        if (isSurfaceMode())
        {
            Info<< "    Surface: " << zoneName
                << ", area: " << areas[zonei] << nl;
        }
        else
        {
            Info<< "    Zone: " << zoneName
                << ", area: " << areas[zonei] << nl;
        }
    }
    Info<< endl;

    if (writeToFile())
    {
        filePtrs_.setSize(zoneNames_.size());

        forAll(filePtrs_, zonei)
        {
            const word& zoneName = zoneNames_[zonei];
            filePtrs_.set(zonei, createFile(zoneName));
            writeFileHeader
            (
                zoneName,
                areas[zonei],
                refDir[zonei],
                filePtrs_[zonei]
            );
        }
    }

    return !zoneNames_.empty();
}


void Foam::functionObjects::fluxSummary::writeFileHeader
(
    const word& zoneName,
    const scalar area,
    const vector& refDir,
    Ostream& os
) const
{
    writeHeader(os, "Flux summary");
    if (isSurfaceMode())
    {
        writeHeaderValue(os, "Surface", zoneName);
    }
    else
    {
        writeHeaderValue(os, "Face zone", zoneName);
    }
    writeHeaderValue(os, "Total area", area);

    switch (mode_)
    {
        case mdFaceZoneAndDirection:
        case mdCellZoneAndDirection:
        case mdSurfaceAndDirection:
        {
            writeHeaderValue(os, "Reference direction", refDir);
            break;
        }
        default:
        {}
    }

    writeHeaderValue(os, "Scale factor", scaleFactor_);

    writeCommented(os, "Time");
    os  << tab << "positive"
        << tab << "negative"
        << tab << "net"
        << tab << "absolute"
        << endl;
}


bool Foam::functionObjects::fluxSummary::execute()
{
    return true;
}


bool Foam::functionObjects::fluxSummary::write()
{
    if (isSurfaceMode())
    {
        return surfaceModeWrite();
    }

    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);

    Log << type() << ' ' << name() << ' '
        << checkFlowType(phi.dimensions(), phi.name()) << " write:" << nl;

    forAll(zoneNames_, zonei)
    {
        const labelList& faceID = faceID_[zonei];
        const labelList& facePatchID = facePatchID_[zonei];
        const boolList&  faceFlips = faceFlip_[zonei];

        scalar phiPos = scalar(0);
        scalar phiNeg = scalar(0);
        scalar phif = scalar(0);

        forAll(faceID, i)
        {
            label facei = faceID[i];
            label patchi = facePatchID[i];

            if (patchi != -1)
            {
                phif = phi.boundaryField()[patchi][facei];
            }
            else
            {
                phif = phi[facei];
            }

            if (faceFlips[i])
            {
                phif *= -1;
            }

            if (phif > 0)
            {
                phiPos += phif;
            }
            else
            {
                phiNeg += phif;
            }
        }

        reduce(phiPos, sumOp<scalar>());
        reduce(phiNeg, sumOp<scalar>());

        phiPos *= scaleFactor_;
        phiNeg *= scaleFactor_;

        scalar netFlux = phiPos + phiNeg;
        scalar absoluteFlux = phiPos - phiNeg;

        Log << "    faceZone " << zoneNames_[zonei] << ':' << nl
            << "        positive : " << phiPos << nl
            << "        negative : " << phiNeg << nl
            << "        net      : " << netFlux << nl
            << "        absolute : " << absoluteFlux
            << nl << endl;

        if (writeToFile())
        {
            filePtrs_[zonei]
                << time_.value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
