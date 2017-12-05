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

#include "polyPDBreader.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyPDBreader, 0);

addToRunTimeSelectionTable(polyConfiguration, polyPDBreader, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPDBreader::polyPDBreader
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
//     const word& name
)
:
    polyConfiguration(molCloud, dict),
    tethered_(false),
    frozen_(false),
    molIds_()
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPDBreader::~polyPDBreader()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyPDBreader::setInitialConfiguration()
{
    // collect bound box information
    vector vMax = vector::zero;
    vector vMin = vector(GREAT, GREAT, GREAT);    
    
    // read in ids
    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        mdInitialiseDict_
    );

    molIds_ = ids.molIds();
    
    const List<word>& idList(molCloud_.cP().molIds());
    
    molIdNames_.setSize(molIds_.size());
    
    forAll(molIdNames_, i)
    {
        molIdNames_[i] = idList[molIds_[i]];
    }
    
    Info << " polyPDBreader on : " << molIdNames_ << endl;
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word regionName(mdInitialiseDict_.lookup("zoneName"));
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("polyPDBreader::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"mdInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    label initialSize = molCloud_.size();
    label nMolsInFile = 0;
    
    if (zone.size())
    {
        Info << "Read positions in zone: " << regionName << endl;
        
        temperature_ = readScalar(mdInitialiseDict_.lookup("temperature"));
        bulkVelocity_ = mdInitialiseDict_.lookup("bulkVelocity");
        
        if (mdInitialiseDict_.found("frozen"))
        {
            frozen_ = Switch(mdInitialiseDict_.lookup("frozen"));
        }        
        
        startPoint_ = mdInitialiseDict_.lookup("displacement");
        
        word name = mdInitialiseDict_.lookup("pdbFileName");  

        ifstream pdbDict(name.c_str());

        if (pdbDict.is_open())
        {

            const reducedUnits& rU = molCloud_.redUnits();
            point globalPosition;
        
            while(pdbDict.good())
            {
                string line;
                getline (pdbDict, line);
                word atomName = line.substr(0,4);
        
                if(atomName == "ATOM")
                {
                    nMolsInFile ++;
                    
//                     string xValue = /*line.substr(30,8)*/;
                    scalar xValue = ::atof((line.substr(30,8)).c_str());
                    globalPosition.x() = xValue/(rU.refLength()*1e10);
                    
//                     word2 = /*line.substr(39,7);*/
                    scalar yValue = ::atof((line.substr(39,7)).c_str());
                    globalPosition.y() = yValue/(rU.refLength()*1e10);
                    
//                     word2 = line.substr(46,8);
                    scalar zValue = ::atof((line.substr(46,8)).c_str());
                    globalPosition.z() = zValue/(rU.refLength()*1e10);
                    
                    
                    testForBoundBox(vMin, vMax, globalPosition);
                    
                    globalPosition += startPoint_;                   
                    
                    // find ID 
                                       
                    label molId = -1;
                    
                   
                    string idNameOneLetterA = line.substr(76,1);
                    string idNameOneLetterB = line.substr(77,1); 
                   
                    string idStr = idNameOneLetterA+idNameOneLetterB;

                    // erase empty spaces
                    idStr.erase(remove_if(idStr.begin(), idStr.end(), isspace), idStr.end());
                    word idName = idStr;
                    label iD = findIndex(molIdNames_, idName);
                    
                    if(iD != -1)
                    {
                        molId = molIds_[iD];

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
                            if (findIndex(zone, cell) != -1)
                            {
                                insertMolecule
                                (
                                    globalPosition,
                                    cell,
                                    tetFace,
                                    tetPt,
                                    molId,
                                    tethered_,
                                    frozen_,
                                    temperature_,
                                    bulkVelocity_
                                );
                            }
                        }
                        else
                        {
                            Info << "position outside = " << globalPosition << endl;
                        }
                    }
                }
            }

            pdbDict.close();
        }
        else
        {
             Info << "Unable to open file"; 
        }
    }

    label finalSize = molCloud_.size();

    nMolsAdded_ = finalSize - initialSize;

    if (Pstream::parRun())
    {
        reduce(nMolsAdded_, sumOp<label>());
    }

    Info << tab << " molecules added: " << nMolsAdded_ << endl;
    
    Info<< nl << tab << " no of molecules in file: " << nMolsInFile 
        << ", molecules not added = " <<  nMolsInFile - nMolsAdded_
        << endl;
    
    
    Info<< nl << "Information on slab of molecules: Vmin = " 
        << vMin << ", vMax = " << vMax 
        << nl << " span: " << vMax - vMin
        <<  endl;
}

        
void polyPDBreader::testForBoundBox
(
    vector& vMin,
    vector& vMax,
    const vector& p
)
{
    if(p.x() > vMax.x())
    {
        vMax.x() = p.x();
    }
    if(p.y() > vMax.y())
    {
        vMax.y() = p.y();
    }    
    if(p.z() > vMax.z())
    {
        vMax.z() = p.z();
    }        
    
    if(p.x() < vMin.x())
    {
        vMin.x() = p.x();
    }
    if(p.y() < vMin.y())
    {
        vMin.y() = p.y();
    }    
    if(p.z() < vMin.z())
    {
        vMin.z() = p.z();
    }  
}




} // End namespace Foam

// ************************************************************************* //
