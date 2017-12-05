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

#include "polyXMOL.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyXMOL, 0);
addToRunTimeSelectionTable(polyField, polyXMOL, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyXMOL::polyXMOL
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    molIds_(),
    excludeSites_(),
    fieldName_(propsDict_.lookup("fieldName")),
    iteration_(0),
    zone_(false),
    regionName_(),
    regionId_(-1)
{
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    if(propsDict_.found("zoneName"))
    {
        zone_ = true;

        const word regionName = propsDict_.lookup("zoneName");
        regionName_ = regionName;

        const cellZoneMesh& cellZones = mesh_.cellZones();

        regionId_ = cellZones.findZoneID(regionName_);

        if(regionId_ == -1)
        {
            FatalErrorIn("polyXMOL::polyXMOL()")
                << "Cannot find region: " << regionName_ << nl << "in: "
                << time_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyXMOL::~polyXMOL()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyXMOL::createField()
{
    selectSiteIds sites
    (
        molCloud_.cP(),
        propsDict_,
        "sitesToExclude"
    );

    List<word> siteNames = sites.siteIdNames();

    excludeSites_.transfer(siteNames);

    Info   << "sites to exclude: " << excludeSites_ << endl;
}

void polyXMOL::calculateField()
{}

void polyXMOL::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        iteration_++;

        write();
    }
}

void polyXMOL::writeInMesh(List<labelField>& molIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        DynamicList<vector> sitePositions(0);
        DynamicList<label> moleculeIds(0);

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                moleculeIds.append(mol().id());

//                 const polyMolecule::constantProperties cP(molCloud_.constProps(mol().id()));

                forAll(mol().sitePositions(), i)
                {
                    if(findIndex(excludeSites_, molCloud_.cP().siteNames(mol().id())[i]) == -1)
                    {
                        sitePositions.append(mol().sitePositions()[i]);
                    }
                }
            }
        }

        //sites[myProc].transfer(sitePositions.shrink());
        //molIds[myProc].transfer(moleculeIds.shrink());

        sites[myProc].transfer(sitePositions);
        molIds[myProc].transfer(moleculeIds);
    }
}

void polyXMOL::writeInZone(List<labelField>& molIds, List<vectorField>& sites)
{
    label myProc =  Pstream::myProcNo();

    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];

    DynamicList<vector> sitePositions(0);
    DynamicList<label> moleculeIds(0);

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                moleculeIds.append(molI->id());

//                 const polyMolecule::constantProperties cP(molCloud_.constProps(molI->id()));

                forAll(molI->sitePositions(), i)
                {
                    if(findIndex(excludeSites_, molCloud_.cP().siteNames(molI->id())[i]) == -1)
                    {
                        sitePositions.append(molI->sitePositions()[i]);
                    }
                }
            }
        }
    }

    //sites[myProc].transfer(sitePositions.shrink());
    //molIds[myProc].transfer(moleculeIds.shrink());

    sites[myProc].transfer(sitePositions);
    molIds[myProc].transfer(moleculeIds);
}

void polyXMOL::write()
{
    List<labelField> molIds(Pstream::nProcs());
    List<vectorField> sites(Pstream::nProcs());

    label myProc =  Pstream::myProcNo();

    if(zone_)
    {
        writeInZone(molIds, sites);
    }
    else
    {
        writeInMesh(molIds, sites);
    }

    label totalSites = sites[myProc].size();
    label totalMols = molIds[myProc].size();

    if (Pstream::parRun())
    {
        reduce(totalSites, sumOp<label>());
        reduce(totalMols, sumOp<label>());
    }

    if (Pstream::parRun())
    {
        // send to master (master does not send)
        if(!Pstream::master())
        {
            const int proc = 0;
            {
                OPstream toNeighbour(Pstream::blocking, proc);
                toNeighbour << sites[myProc] << molIds[myProc];
            }
        }

        //- receiving (master only receives)
        if(Pstream::master())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    vectorField sitesProc;
                    labelField molIdsProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> sitesProc >> molIdsProc;
                    }

                    sites[p].setSize(sitesProc.size());
                    molIds[p].setSize(molIdsProc.size());

                    forAll(sitesProc, i)
                    {
                        sites[p][i] = sitesProc[i];
                    }

                    forAll(molIdsProc, i)
                    {
                        molIds[p][i] = molIdsProc[i];
                    }
                }
            }
        }
    }

    if(Pstream::master())
    {
        const reducedUnits& rU = molCloud_.redUnits();

        fileName fName(timePath_/"polyMoleculeCloud_"+fieldName_+".xmol");
    
        OFstream str(fName);
    
        str << totalSites << nl << "polyMoleculeCloud site positions in angstroms" << nl;

        // for all processors
        forAll(molIds, p)
        {
            label posCounter = 0;

            forAll(molIds[p], i)
            {
                label molId = molIds[p][i];

//                 const polyMolecule::constantProperties cP(molCloud_.constProps(molId));

                label n = molCloud_.cP().nSites(molId);
            
                for (int i = 0; i < n; i++)
                {
                    if(findIndex(excludeSites_, molCloud_.cP().siteNames(molId)[i]) == -1)
                    {
                        vector rS = sites[p][posCounter]*rU.refLength()*1e10;

                        str <<  molCloud_.cP().siteNames(molId)[i]
                            << ' ' << rS.x()
                            << ' ' << rS.y()
                            << ' ' << rS.z()
                            << nl;

                        posCounter++;
                    }
                }
            }
        }
    }
}

void polyXMOL::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyXMOL::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyXMOL::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
