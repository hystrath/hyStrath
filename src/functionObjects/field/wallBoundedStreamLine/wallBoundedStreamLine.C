/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "wallBoundedStreamLine.H"
#include "wallBoundedStreamLineParticleCloud.H"
#include "sampledSet.H"
#include "faceSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallBoundedStreamLine, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallBoundedStreamLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tetIndices Foam::functionObjects::wallBoundedStreamLine::findNearestTet
(
    const PackedBoolList& isWallPatch,
    const point& seedPt,
    const label celli
) const
{
    const cell& cFaces = mesh_.cells()[celli];

    label minFacei = -1;
    label minTetPti = -1;
    scalar minDistSqr = sqr(GREAT);

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        if (isWallPatch[facei])
        {
            const face& f = mesh_.faces()[facei];
            const label fp0 = mesh_.tetBasePtIs()[facei];
            const point& basePoint = mesh_.points()[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = mesh_.points()[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = mesh_.points()[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);

                scalar d2 = magSqr(tri.centre() - seedPt);
                if (d2 < minDistSqr)
                {
                    minDistSqr = d2;
                    minFacei = facei;
                    minTetPti = i-1;
                }
                fp = nextFp;
            }
        }
    }

    // Put particle in tet
    return tetIndices
    (
        celli,
        minFacei,
        minTetPti,
        mesh_
    );
}


void Foam::functionObjects::wallBoundedStreamLine::track()
{
    // Determine the 'wall' patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are the faces that need to be followed

    autoPtr<indirectPrimitivePatch> boundaryPatch(wallPatch());
    PackedBoolList isWallPatch(mesh_.nFaces());
    forAll(boundaryPatch().addressing(), i)
    {
        isWallPatch[boundaryPatch().addressing()[i]] = 1;
    }



    // Find nearest wall particle for the seedPoints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IDLList<wallBoundedStreamLineParticle> initialParticles;
    wallBoundedStreamLineParticleCloud particles
    (
        mesh_,
        cloudName_,
        initialParticles
    );

    {
        // Get the seed points
        // ~~~~~~~~~~~~~~~~~~~

        const sampledSet& seedPoints = sampledSetPtr_();


        forAll(seedPoints, i)
        {
            const point& seedPt = seedPoints[i];
            label celli = seedPoints.cells()[i];

            tetIndices ids(findNearestTet(isWallPatch, seedPt, celli));

            if (ids.face() != -1 && isWallPatch[ids.face()])
            {
                //Pout<< "Seeding particle :" << nl
                //    << "     seedPt:" << seedPt << nl
                //    << "     face  :" << ids.face() << nl
                //    << "     at    :" << mesh_.faceCentres()[ids.face()] << nl
                //    << "     cell  :" << mesh_.cellCentres()[ids.cell()] << nl
                //    << endl;

                particles.addParticle
                (
                    new wallBoundedStreamLineParticle
                    (
                        mesh_,
                        ids.faceTri(mesh_).centre(),
                        ids.cell(),
                        ids.face(),     // tetFace
                        ids.tetPt(),
                        -1,             // not on a mesh edge
                        -1,             // not on a diagonal edge
                        lifeTime_       // lifetime
                    )
                );
            }
            else
            {
                Pout<< type() << " : ignoring seed " << seedPt
                    << " since not in wall cell." << endl;
            }
        }
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    Log << type() << " : seeded " << nSeeds << " particles." << endl;



    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar>> vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector>> vvInterp;

    label UIndex = -1;

    initInterpolations
    (
        nSeeds,
        UIndex,
        vsFlds,
        vsInterp,
        vvFlds,
        vvInterp
    );

    // Additional particle info
    wallBoundedStreamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        trackLength_,   // fixed track length
        isWallPatch,    // which faces are to follow

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoundedStreamLine::wallBoundedStreamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    streamLineBase(name, runTime, dict)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoundedStreamLine::~wallBoundedStreamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallBoundedStreamLine::read(const dictionary& dict)
{
    if (streamLineBase::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;


        // Make sure that the mesh is trackable
        if (debug)
        {
            // 1. Positive volume decomposition tets
            faceSet faces(mesh_, "lowQualityTetFaces", mesh_.nFaces()/100+1);
            if
            (
                polyMeshTetDecomposition::checkFaceTets
                (
                    mesh_,
                    polyMeshTetDecomposition::minTetQuality,
                    true,
                    &faces
                )
            )
            {
                label nFaces = returnReduce(faces.size(), sumOp<label>());

                WarningInFunction
                    << "Found " << nFaces
                    <<" faces with low quality or negative volume "
                    << "decomposition tets. Writing to faceSet " << faces.name()
                    << endl;
            }

            // 2. All edges on a cell having two faces
            EdgeMap<label> numFacesPerEdge;
            forAll(mesh_.cells(), celli)
            {
                const cell& cFaces = mesh_.cells()[celli];

                numFacesPerEdge.clear();

                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    const face& f = mesh_.faces()[facei];
                    forAll(f, fp)
                    {
                        const edge e(f[fp], f.nextLabel(fp));
                        EdgeMap<label>::iterator eFnd =
                            numFacesPerEdge.find(e);
                        if (eFnd != numFacesPerEdge.end())
                        {
                            eFnd()++;
                        }
                        else
                        {
                            numFacesPerEdge.insert(e, 1);
                        }
                    }
                }

                forAllConstIter(EdgeMap<label>, numFacesPerEdge, iter)
                {
                    if (iter() != 2)
                    {
                        FatalErrorInFunction
                            << "problem cell:" << celli
                            << abort(FatalError);
                    }
                }
            }
        }
    }

    return true;
}


// ************************************************************************* //
