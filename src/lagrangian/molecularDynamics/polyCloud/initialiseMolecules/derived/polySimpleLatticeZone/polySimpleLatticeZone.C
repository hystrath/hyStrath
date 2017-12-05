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

#include "polySimpleLatticeZone.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polySimpleLatticeZone, 0);
addToRunTimeSelectionTable(polyConfiguration, polySimpleLatticeZone, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polySimpleLatticeZone::polySimpleLatticeZone
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyConfiguration(molCloud, dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polySimpleLatticeZone::~polySimpleLatticeZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polySimpleLatticeZone::setInitialConfiguration()
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word regionName(mdInitialiseDict_.lookup("zoneName"));
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("polySimpleLatticeZone::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"mdInitialiseDict"
            << exit(FatalError);
    }

    const labelList& fillZone = cellZones[zoneId];

    label initialSize = molCloud_.size();

    Info << nl << "Quick lattice in zone: " << regionName << endl;

    const scalar temperature(readScalar(mdInitialiseDict_.lookup("temperature")));

    const vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));

    const scalar minSpacing(readScalar(mdInitialiseDict_.lookup("minSpacing")));

    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polySimpleLatticeZone::setInitialConfiguration()")
            << "Cannot find molecule id: " << molIdName << nl << "in idList."
            << exit(FatalError);
    }

    bool frozen = false;

    if (mdInitialiseDict_.found("frozen"))
    {
        frozen = Switch(mdInitialiseDict_.lookup("frozen"));
    }

    bool tethered = false;

    if (mdInitialiseDict_.found("tethered"))
    {
        tethered = Switch(mdInitialiseDict_.lookup("tethered"));
    }

    scalar offset = 0.0;

    if (mdInitialiseDict_.found("offset"))
    {
        offset = readScalar
        (
            mdInitialiseDict_.lookup("offset")
        );
    }

    vector offsetMin = vector::zero;

    if (mdInitialiseDict_.found("offsetMin"))
    {
        offsetMin = (mdInitialiseDict_.lookup("offsetMin"));
    }

    vector offsetMax = vector::zero;

    if (mdInitialiseDict_.found("offsetMax"))
    {
        offsetMax = (mdInitialiseDict_.lookup("offsetMax"));
    }

    vector displacement = vector::zero;

    if (mdInitialiseDict_.found("displacement"))
    {
        displacement =
        (
            mdInitialiseDict_.lookup("displacement")
        );
    }

    bool displaceOutOfZone = false;

    if (mdInitialiseDict_.found("displaceOutOfZone"))
    {
        displaceOutOfZone = Switch(mdInitialiseDict_.lookup("displaceOutOfZone"));
    }

    bool distributeSpacing = false;

    if (mdInitialiseDict_.found("distributeSpacing"))
    {
        distributeSpacing = Switch(mdInitialiseDict_.lookup("distributeSpacing"));
    }

    scalar zoneVolume = 0.0;

    forAll(fillZone, c)
    {
        const label& cellI = fillZone[c];
        zoneVolume += mesh_.cellVolumes()[cellI];
    }

    if(Pstream::parRun())
    {
        reduce(zoneVolume, sumOp<scalar>());
    }

    if(zoneVolume == 0.0)
    {
        FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
            << "zone volume is zero: "
            << regionName
            << abort(FatalError);
    }

    Info << "volume of zone: " << zoneVolume << endl;

    vector zoneMin = VGREAT*vector::one;

    vector zoneMax = -VGREAT*vector::one;

    forAll (fillZone, c)
    {
        const label& cellI = fillZone[c];

        const labelList& cellPoints = mesh_.cellPoints()[cellI];

        forAll(cellPoints, cP)
        {
            const label& p = cellPoints[cP];

            const point& pointI = mesh_.points()[p];

            if (pointI.x() > zoneMax.x())
            {
                zoneMax.x() = pointI.x();
            }
            if (pointI.x() < zoneMin.x())
            {
                zoneMin.x() = pointI.x();
            }
            if (pointI.y() > zoneMax.y())
            {
                zoneMax.y() = pointI.y();
            }
            if (pointI.y() < zoneMin.y())
            {
                zoneMin.y() = pointI.y();
            }
            if (pointI.z() > zoneMax.z())
            {
                zoneMax.z() = pointI.z();
            }
            if (pointI.z() < zoneMin.z())
            {
                zoneMin.z() = pointI.z();
            }
        }
    }

    if (Pstream::parRun())
    {
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << zoneMin << zoneMax ;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                vector zoneMinProc;
                vector zoneMaxProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> zoneMinProc >> zoneMaxProc;
                }

                if (zoneMaxProc.x() > zoneMax.x())
                {
                    zoneMax.x() = zoneMaxProc.x();
                }
                if (zoneMaxProc.y() > zoneMax.y())
                {
                    zoneMax.y() = zoneMaxProc.y();
                }
                if (zoneMaxProc.z() > zoneMax.z())
                {
                    zoneMax.z() = zoneMaxProc.z();
                }


                if (zoneMinProc.x() < zoneMin.x())
                {
                    zoneMin.x() = zoneMinProc.x();
                }
                if (zoneMinProc.y() < zoneMin.y())
                {
                    zoneMin.y() = zoneMinProc.y();
                }
                if (zoneMinProc.z() < zoneMin.z())
                {
                    zoneMin.z() = zoneMinProc.z();
                }
            }
        }
    }

    if(mag(offset) > 0.0)
    {
        Info << "Original min: " << zoneMin << " max: " << zoneMax << endl;

        vector unitV = (zoneMax - zoneMin)/mag(zoneMax-zoneMin);
        zoneMin -= unitV*offset;
        zoneMax += unitV*offset;
    }

    if(mag(offsetMin) > 0.0)
    {
       zoneMin += offsetMin;
    }

    if(mag(offsetMax) > 0.0)
    {
       zoneMax += offsetMax;
    }

    Info << "min: " << zoneMin << ", max: " << zoneMax << endl;

    boundBox bb(zoneMin, zoneMax);

    scalar bbVol = bb.span().x() * bb.span().y() * bb.span().z();

    Info << "volume of bound box: " << bbVol << endl;

    scalar numberDensity = 0.0;

    if (mdInitialiseDict_.found("numberDensity"))
    {
        scalar rho = readScalar
        (
            mdInitialiseDict_.lookup("numberDensity")
        );

        numberDensity = rho;

        if (numberDensity < VSMALL)
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "numberDensity too small, not filling zone "
                << regionName
                << abort(FatalError);
        }
    }
    else if (mdInitialiseDict_.found("massDensity"))
    {
//         const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));

        scalar mass = molCloud_.cP().mass(molId);

        Info << "mass: " << mass << endl;

        scalar massDensity = readScalar
        (
            mdInitialiseDict_.lookup("massDensity")
        );

        numberDensity = massDensity / mass;

        if (massDensity < VSMALL)
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity too small, not filling zone "
                << regionName
                << abort(FatalError);
        }
    }
    else if (mdInitialiseDict_.found("massDensitySI"))
    {
//         const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));

        scalar mass = molCloud_.cP().mass(molId);

        Info << "mass: " << mass << endl;

        scalar massDensity = readScalar
        (
            mdInitialiseDict_.lookup("massDensitySI")
        );
        scalar m = molCloud_.redUnits().refMass();
        scalar l = molCloud_.redUnits().refLength();
        
        numberDensity = massDensity*l*l*l / (m*mass);

        if (massDensity < VSMALL)
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity too small, not filling zone "
                << regionName
                << abort(FatalError);
        }
    }    
    else
    {
        FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
            << "massDensity or numberDensity not specified " << nl
            << abort(FatalError);
    }

    scalar nMols = numberDensity * bbVol;

    Info << "number of molecules to insert (scalar): " << nMols << endl;

    scalar spacing = pow(  (1.0/numberDensity), (1.0/3.0) );

    Info << "spacing: " << spacing << endl;

    if(spacing < minSpacing)
    {
        spacing = minSpacing;
    }

    scalar spacingX = spacing;
    scalar spacingY = spacing;
    scalar spacingZ = spacing;

    label nMolsX = label( bb.span().x()/spacing);
    label nMolsY = label( bb.span().y()/spacing);
    label nMolsZ = label( bb.span().z()/spacing);

    scalar spacingResX = bb.span().x()-(nMolsX*spacing);
    scalar spacingResY = bb.span().y()-(nMolsY*spacing);
    scalar spacingResZ = bb.span().z()-(nMolsZ*spacing);

    if(distributeSpacing)
    {
        spacingResX = 0.0;
        spacingResY = 0.0;
        spacingResZ = 0.0;

        spacingX = (bb.span().x()/nMolsX);
        spacingY = (bb.span().y()/nMolsY);
        spacingZ = (bb.span().z()/nMolsZ);

        Info<< "distributing spacing to: " 
            << " spacingX: " << spacingX
            << ", spacingY: " << spacingY 
            << ", spacingZ: " << spacingZ
            << endl;
    }

    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;

    Info << "nMolsX: " << nMolsX << ", nMolsY: " <<  nMolsY 
         << ", nMolsZ: " <<  nMolsZ << endl; 

    vector globalPosition = vector::zero;

    for (label i = 0; i < nMolsX; i++)
    {
        x = (0.5 + i)*spacingX + spacingResX*0.5;

        for (label j = 0; j < nMolsY; j++)
        {
            y = (0.5 + j)*spacingY + spacingResY*0.5;

            for (label k = 0; k < nMolsZ; k++)
            {
                z = (0.5 + k)*spacingZ + spacingResZ*0.5;

                globalPosition = bb.min() + vector(x, y, z) + displacement;
            
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

                if(cell != -1)
                {
                    if( (findIndex(fillZone, cell) != -1) || (displaceOutOfZone) )
                    {
                        insertMolecule
                        (
                            globalPosition,
                            cell,
                            tetFace,
                            tetPt,
                            molId,
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
