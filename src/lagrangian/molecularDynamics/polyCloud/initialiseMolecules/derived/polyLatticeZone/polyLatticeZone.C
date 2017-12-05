/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "polyLatticeZone.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyLatticeZone, 0);

addToRunTimeSelectionTable(polyConfiguration, polyLatticeZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyLatticeZone::polyLatticeZone
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict/*, name*/)
//     propsDict_(dict.subDict(typeName + "Properties"))
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyLatticeZone::~polyLatticeZone()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyLatticeZone::setInitialConfiguration()
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word regionName(mdInitialiseDict_.lookup("zoneName"));
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("atomisticLatticeZone::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"mdInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    label initialSize = molCloud_.size();


    if (zone.size())
    {
        Info << "Lattice in zone: " << regionName << endl;

        const scalar temperature
        (
            readScalar(mdInitialiseDict_.lookup("temperature"))
        );

        const vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));

        List<word> latticeIds
        (
            mdInitialiseDict_.lookup("latticeIds")
        );


        const List<word>& idList(molCloud_.cP().molIds());

        // test for lattice ids in idList
        forAll(latticeIds, i)
        {
            label molId = findIndex(idList, latticeIds[i]);

            if(molId == -1)
            {
                FatalErrorIn("polyLatticeZone::setInitialConfiguration()")
                    << "Cannot find molecule id: " << latticeIds[i] << nl << "in idList."
                    << exit(FatalError);
            }
        }

        List<vector> latticePositions
        (
            mdInitialiseDict_.lookup("latticePositions")
        );

        if(latticeIds.size() != latticePositions.size())
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "latticeIds and latticePositions must be the same "
                << " size." << nl
                << abort(FatalError);
        }

        diagTensor latticeCellShape
        (
            mdInitialiseDict_.lookup("latticeCellShape")
        );

        scalar latticeCellScale = 0.0;

        if (mdInitialiseDict_.found("numberDensity"))
        {
            scalar numberDensity = readScalar
            (
                mdInitialiseDict_.lookup("numberDensity")
            );

            if (numberDensity < VSMALL)
            {
                FatalErrorIn("polyMoleculeCloud::initialiseMolecules")
                    << "numberDensity too small, not filling zone "
                    << zone.name()
                    << abort(FatalError);
            }

            Info << "latticeIds.size(): " << latticeIds.size() << endl;

            latticeCellScale = pow
            (
                latticeIds.size()/(det(latticeCellShape)*numberDensity),
                (1.0/3.0)
            );
        }
        else if (mdInitialiseDict_.found("massDensity"))
        {
            scalar unitCellMass = 0.0;

            forAll(latticeIds, i)
            {
                label id = findIndex(molCloud_.cP().molIds(), latticeIds[i]);

//                 const polyMolecule::constantProperties& cP(molCloud_.constProps(id));

                unitCellMass += molCloud_.cP().mass(id);
            }

            Info << "unitCellMass: " <<unitCellMass << endl;

            scalar massDensity = readScalar
            (
                mdInitialiseDict_.lookup("massDensity")
            );

            if (massDensity < VSMALL)
            {
                FatalErrorIn("polyMoleculeCloud::initialiseMolecules")
                    << "massDensity too small, not filling zone "
                    << zone.name()
                    << abort(FatalError);
            }


            latticeCellScale = pow
            (
                unitCellMass/(det(latticeCellShape)*massDensity),
                (1.0/3.0)
            );
        }
        else
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity or numberDensity not specified " << nl
                << abort(FatalError);
        }

        latticeCellShape *= latticeCellScale;

        vector anchor(mdInitialiseDict_.lookup("anchor"));

        bool tethered = false;

        if (mdInitialiseDict_.found("tetherSiteIds"))
        {
            tethered = bool
            (
                List<word>(mdInitialiseDict_.lookup("tetherSiteIds")).size()
            );
        }

        bool frozen = false;//*****

        if (mdInitialiseDict_.found("frozenSiteIds"))
        {
            frozen = bool
            (
                List<word>(mdInitialiseDict_.lookup("frozenSiteIds")).size()
            );
        }

        const vector orientationAngles
        (
            mdInitialiseDict_.lookup("orientationAngles")
        );

        scalar phi
        (
            orientationAngles.x()*constant::mathematical::pi/180.0
        );

        scalar theta
        (
            orientationAngles.y()*constant::mathematical::pi/180.0
        );

        scalar psi
        (
            orientationAngles.z()*constant::mathematical::pi/180.0
        );

        const tensor R
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );

        // Find the optimal anchor position.  Finding the approximate
        // mid-point of the zone of cells and snapping to the nearest
        // lattice location.

        vector zoneMin = VGREAT*vector::one;

        vector zoneMax = -VGREAT*vector::one;

        forAll (zone, cell)
        {
            const point cellCentre = mesh_.cellCentres()[zone[cell]];

            if (cellCentre.x() > zoneMax.x())
            {
                zoneMax.x() = cellCentre.x();
            }
            if (cellCentre.x() < zoneMin.x())
            {
                zoneMin.x() = cellCentre.x();
            }
            if (cellCentre.y() > zoneMax.y())
            {
                zoneMax.y() = cellCentre.y();
            }
            if (cellCentre.y() < zoneMin.y())
            {
                zoneMin.y() = cellCentre.y();
            }
            if (cellCentre.z() > zoneMax.z())
            {
                zoneMax.z() = cellCentre.z();
            }
            if (cellCentre.z() < zoneMin.z())
            {
                zoneMin.z() = cellCentre.z();
            }
        }

        point zoneMid = 0.5*(zoneMin + zoneMax);

        point latticeMid = inv(latticeCellShape) & (R.T()
        & (zoneMid - anchor));

        point latticeAnchor
        (
            label(latticeMid.x() + 0.5*sign(latticeMid.x())),
            label(latticeMid.y() + 0.5*sign(latticeMid.y())),
            label(latticeMid.z() + 0.5*sign(latticeMid.z()))
        );

        anchor += (R & (latticeCellShape & latticeAnchor));

        // Continue trying to place molecules as long as at
        // least one molecule is placed in each iteration.
        // The "|| totalZoneMols == 0" condition means that the
        // algorithm will continue if the origin is outside the
        // zone.

        label n = 0;

        label totalZoneMols = 0;

        label molsPlacedThisIteration = 0;

        while
        (
            molsPlacedThisIteration != 0
            || totalZoneMols == 0
        )
        {
            label sizeBeforeIteration = molCloud_.size();

            bool partOfLayerInBounds = false;

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // start of placement of molecules
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (n == 0)
            {
                // Special treatment is required for the first position,
                // i.e. iteration zero.

                labelVector unitCellLatticePosition(0,0,0);

                forAll(latticePositions, p)
                {
                    label id = findIndex(molCloud_.cP().molIds(), latticeIds[p]);

                    const vector& latticePosition =
                    vector
                    (
                        unitCellLatticePosition.x(),
                        unitCellLatticePosition.y(),
                        unitCellLatticePosition.z()
                    )
                    + latticePositions[p];

                    point globalPosition =
                    anchor + (R & (latticeCellShape & latticePosition));

                    partOfLayerInBounds =
                    mesh_.bounds().contains(globalPosition);


                    label cell = -1;
                    label tetFace = -1;
                    label tetPt = -1;

                    mesh_.findCellFacePt
                    (
                        globalPosition,
                        cell,
                        tetFace,
                        tetPt
                    );

                    if (findIndex(zone, cell) != -1)
                    {
                        insertMolecule
                        (
                            globalPosition,
                            cell,
                            tetFace,
                            tetPt,                         
                            id,
                            tethered,
                            frozen,
                            temperature,
                            bulkVelocity
                        );
                    }
                }
            }
            else
            {
                // Place top and bottom caps.

                labelVector unitCellLatticePosition(0,0,0);

                for
                (
                    unitCellLatticePosition.z() = -n;
                    unitCellLatticePosition.z() <= n;
                    unitCellLatticePosition.z() += 2*n
                )
                {
                    for
                    (
                        unitCellLatticePosition.y() = -n;
                        unitCellLatticePosition.y() <= n;
                        unitCellLatticePosition.y()++
                    )
                    {
                        for
                        (
                            unitCellLatticePosition.x() = -n;
                            unitCellLatticePosition.x() <= n;
                            unitCellLatticePosition.x()++
                        )
                        {
                            forAll(latticePositions, p)
                            {
                                label id = findIndex
                                (
                                    molCloud_.cP().molIds(),
                                    latticeIds[p]
                                );

                                const vector& latticePosition =
                                vector
                                (
                                    unitCellLatticePosition.x(),
                                    unitCellLatticePosition.y(),
                                    unitCellLatticePosition.z()
                                )
                                + latticePositions[p];

                                point globalPosition = anchor
                                + (R
                                &(latticeCellShape & latticePosition));

                                partOfLayerInBounds =
                                mesh_.bounds().contains(globalPosition);

                                label cell = -1;
                                label tetFace = -1;
                                label tetPt = -1;

                                mesh_.findCellFacePt
                                (
                                    globalPosition,
                                    cell,
                                    tetFace,
                                    tetPt
                                );

                                if (findIndex(zone, cell) != -1)
                                {
                                    insertMolecule
                                    (
                                        globalPosition,
                                        cell,
                                        tetFace,
                                        tetPt,
                                        id,
                                        tethered,
                                        frozen,
                                        temperature,
                                        bulkVelocity
                                    );
                                }
                            }
                        }
                    }
                }

                for
                (
                    unitCellLatticePosition.z() = -(n-1);
                    unitCellLatticePosition.z() <= (n-1);
                    unitCellLatticePosition.z()++
                )
                {
                    for (label iR = 0; iR <= 2*n -1; iR++)
                    {
                        unitCellLatticePosition.x() = n;

                        unitCellLatticePosition.y() = -n + (iR + 1);

                        for (label iK = 0; iK < 4; iK++)
                        {
                            forAll(latticePositions, p)
                            {
                                label id = findIndex
                                (
                                    molCloud_.cP().molIds(),
                                    latticeIds[p]
                                );

                                const vector& latticePosition =
                                vector
                                (
                                    unitCellLatticePosition.x(),
                                    unitCellLatticePosition.y(),
                                    unitCellLatticePosition.z()
                                )
                                + latticePositions[p];

                                point globalPosition = anchor
                                + (R
                                &(latticeCellShape & latticePosition));

                                partOfLayerInBounds =
                                mesh_.bounds().contains(globalPosition);

                                label cell = -1;
                                label tetFace = -1;
                                label tetPt = -1;

                                mesh_.findCellFacePt
                                (
                                    globalPosition,
                                    cell,
                                    tetFace,
                                    tetPt
                                );

                                if (findIndex(zone, cell) != -1)
                                {
                                    insertMolecule
                                    (
                                        globalPosition,
                                        cell,
                                        tetFace,
                                        tetPt,
                                        id,
                                        tethered,
                                        frozen,
                                        temperature,
                                        bulkVelocity
                                    );
                                }
                            }

                            unitCellLatticePosition =
                            labelVector
                            (
                                - unitCellLatticePosition.y(),
                                unitCellLatticePosition.x(),
                                unitCellLatticePosition.z()
                            );
                        }
                    }
                }
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // end of placement of molecules
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if
            (
                totalZoneMols == 0
                && !partOfLayerInBounds
            )
            {
                WarningIn("Foam::polyMoleculeCloud::initialiseMolecules()")
                    << "A whole layer of unit cells was placed "
                    << "outside the bounds of the mesh, but no "
                    << "molecules have been placed in zone '"
                    << zone.name()
                    << "'.  This is likely to be because the zone "
                    << "has few cells ("
                    << zone.size()
                    << " in this case) and no lattice position "
                    << "fell inside them.  "
                    << "Aborting filling this zone."
                    << endl;

                break;
            }

            molsPlacedThisIteration = molCloud_.size() - sizeBeforeIteration;

            totalZoneMols += molsPlacedThisIteration;

            n++;
        }
    }

    label finalSize = molCloud_.size();

    nMolsAdded_ = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nMolsAdded_, sumOp<label>());
    }

    Info << tab << " molecules added: " << nMolsAdded_ << endl;

}



} // End namespace Foam

// ************************************************************************* //
