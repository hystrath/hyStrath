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
    multiZoneBins

Description

\*----------------------------------------------------------------------------*/

#include "multiZoneBins.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(multiZoneBins, 0);

addToRunTimeSelectionTable(binModel, multiZoneBins, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
multiZoneBins::multiZoneBins
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    binModel(mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    startPoint_(propsDict_.lookup("startPoint")),
    endPoint_(propsDict_.lookup("endPoint")),
    unitVector_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    nBins_(0),
    regionIds_(),
    volumes_(),
    binWidths_(),
    cellRegionAddressing_()
{
    // read in list of zone names
    const List<word> zoneNames (propsDict_.lookup("zoneNames"));

    if(!zoneNames.size())
    {
        FatalErrorIn("multiZoneBins::multiZoneBins()")
            << "Define zoneNames " << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    // check if a cell zone is defined more than once in the input list
    {
        DynamicList<word> regionNames(0);

        forAll(zoneNames, i)
        {
            const word& zoneName(zoneNames[i]);
    
            if(findIndex(regionNames, zoneName) == -1)
            {
                regionNames.append(zoneName);
            }
            else
            {
                FatalErrorIn("multiZoneBins::multiZoneBins()")
                    << "Zone name: " << zoneName << " cannot be defined twice "
                    << nl << "in: "
                    << mesh_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);
            }
        }
    }


    // test if specified region names belong to the mesh's cell zones

    const cellZoneMesh& cellZones = mesh_.cellZones();

    nBins_ = zoneNames.size();

    regionIds_.setSize(nBins_, -1);

    forAll(regionIds_, r)
    {
        regionIds_[r] = cellZones.findZoneID(zoneNames[r]);
    
        if(regionIds_[r] == -1)
        {
            FatalErrorIn("multiZoneBins::multiZoneBins()")
                << "Cannot find region: " << zoneNames[r] 
                << " in mesh's cell zones."
                << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }


    //-set the total volumes

    volumes_.setSize(nBins_, 0.0);

    forAll(regionIds_, r)
    {
        const labelList& cells = cellZones[regionIds_[r]];

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            volumes_[r] += mesh_.cellVolumes()[cellI];
        }

        if (Pstream::parRun())
        {
            reduce(volumes_[r], sumOp<scalar>());
        }
    }

    // set cell-region adressing
    // check for overlapping cells too

    cellRegionAddressing_.setSize(mesh_.nCells(), -1);

    forAll(regionIds_, r)
    {
        const labelList& cells = cellZones[regionIds_[r]];

        forAll(cells, c)
        {
            const label& cellI = cells[c];

            if(cellRegionAddressing_[cellI] == -1)
            {
                cellRegionAddressing_[cellI] = r;
            }
            else
            {
                FatalErrorIn("multiZoneBins::multiZoneBins()")
                    << "There is an overlap present in cell-zone: " << zoneNames[r] 
                    << " with zone: " << zoneNames[cellRegionAddressing_[cellI]]
                    << nl << "in: "
                    << mesh_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);
            }
        }
    }

    // set the binWidths

    binWidths_.setSize(nBins_, 0.0);

    scalar rSEMag = mag(endPoint_ - startPoint_);

    forAll(regionIds_, r)
    {
        const labelList& cells = cellZones[regionIds_[r]];

        scalar rMax = 0.0;
        scalar rMin = GREAT;

        forAll(cells, c)
        {
            const label& cellI = cells[c];

            const labelList& points = mesh_.cellPoints()[cellI];

            forAll(points, p)
            {
                vector rI = mesh_.points()[points[p]];

                scalar nD = (rI - startPoint_) & unitVector_;

                if((nD >= 0) && (nD <= rSEMag))
                {
                    if(nD > rMax)
                    {
                        rMax = nD;
                    }
                    if(nD < rMin)
                    {
                        rMin = nD;
                    }
                }
            }
        }

        //- parallel communication
        if(Pstream::parRun())
        {
            //-sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << rMax << rMin;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar rMaxProc;
                    scalar rMinProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> rMaxProc >> rMinProc;
                    }
        
                    if(rMaxProc > rMax)
                    {
                        rMax = rMaxProc;
                    }

                    if(rMinProc < rMin)
                    {
                        rMin = rMinProc;
                    }
                }
            }
        }

        binWidths_[r] = rMax - rMin;

        Info<<"rMin : " << rMin << ", rMax: " << rMax
            << " binWidth: " << binWidths_[r]
            << endl;
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

multiZoneBins::~multiZoneBins()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


label multiZoneBins::isPointWithinBin
(
    const vector& rI,
    const label& cellI
)
{
    label binNumber = cellRegionAddressing_[cellI];

    return binNumber;
}


scalarField multiZoneBins::binPositions()
{
    scalarField positions(nBins_, 0.0);

    if(nBins_ > 0)
    {
        positions[0] = 0.5*binWidths_[0];
    
        for (label i = 1; i < nBins_; i++)
        {
            positions[i] = positions[i-1] + 0.5*binWidths_[i-1] + 0.5*binWidths_[i];
        }
    }

//     Info << "positions: " << positions << endl;

    return positions;
}

vectorField multiZoneBins::bins()
{
    vectorField positions(nBins_, vector::zero);

    if(nBins_ > 0)
    {
        positions[0] = startPoint_ + 0.5*binWidths_[0]*unitVector_;

        for (label i = 1; i < nBins_; i++)
        {
            positions[i] = positions[i-1] + 0.5*binWidths_[i-1]*unitVector_ 
                                          + 0.5*binWidths_[i]*unitVector_;
        }
    }

//     Info << "positions: " << positions << endl;

    return positions;
}

const label& multiZoneBins::nBins() const
{
    return nBins_;
}

scalar multiZoneBins::binVolume(const label& n)
{
    return volumes_[n];
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
