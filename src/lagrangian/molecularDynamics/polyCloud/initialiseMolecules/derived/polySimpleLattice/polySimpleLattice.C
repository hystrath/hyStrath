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

#include "polySimpleLattice.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polySimpleLattice, 0);

addToRunTimeSelectionTable(polyConfiguration, polySimpleLattice, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polySimpleLattice::polySimpleLattice
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

polySimpleLattice::~polySimpleLattice()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polySimpleLattice::setInitialConfiguration()
{
    label initialSize = molCloud_.size();

    Info << nl << "Building quick simple lattice: " << endl;

    const scalar temperature(readScalar(mdInitialiseDict_.lookup("temperature")));

    const vector bulkVelocity(mdInitialiseDict_.lookup("bulkVelocity"));

    bool computeSpacingFromDensity = false;

    if (mdInitialiseDict_.found("computeSpacingFromDensity"))
    {
        computeSpacingFromDensity = Switch(mdInitialiseDict_.lookup("computeSpacingFromDensity"));
    }

    scalar spacing = 0.0;

    if(!computeSpacingFromDensity)
    {
        spacing = readScalar(mdInitialiseDict_.lookup("spacing"));
    }

    const word molIdName(mdInitialiseDict_.lookup("molId")); 
    const List<word>& idList(molCloud_.cP().molIds());

    label molId = findIndex(idList, molIdName);

    if(molId == -1)
    {
        FatalErrorIn("polySimpleLattice::setInitialConfiguration()")
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

    bool deleteOutOfMesh = true;

    if (mdInitialiseDict_.found("deleteOutOfMesh"))
    {
        deleteOutOfMesh = Switch(mdInitialiseDict_.lookup("deleteOutOfMesh"));
    }

    bool deleteOutOfZone = false;
    label zoneId = -1;
    labelList cells;

    if (mdInitialiseDict_.found("deleteOutOfZone"))
    {
        deleteOutOfZone = Switch(mdInitialiseDict_.lookup("deleteOutOfZone"));

        if(deleteOutOfZone)
        {
            const cellZoneMesh& cellZones = mesh_.cellZones();
            const word regionName(mdInitialiseDict_.lookup("zoneName"));
            zoneId = cellZones.findZoneID(regionName);

            if(zoneId == -1)
            {
                FatalErrorIn("polySimpleLatticeZone::setInitialConfiguration()")
                    << "Cannot find region: " << regionName << nl << "in: "
                    << mesh_.time().system()/"mdInitialiseDict"
                    << exit(FatalError);
            }

            const labelList& fillZone = cellZones[zoneId];
            cells.setSize(fillZone.size());

            forAll(fillZone, c)
            {
                cells[c] = fillZone[c];
            }
        }
    }

    bool adjustSpacing = true;

    if (mdInitialiseDict_.found("adjustSpacing"))
    {
        adjustSpacing = Switch(mdInitialiseDict_.lookup("adjustSpacing"));
    }

    bool deleteOverlaps = false;
    
    scalar rOv = spacing;

    if (mdInitialiseDict_.found("deleteOverlaps"))
    {
        deleteOverlaps = Switch(mdInitialiseDict_.lookup("deleteOverlaps"));
        rOv = readScalar(mdInitialiseDict_.lookup("overlapSpacing"));
    }

    const vector startPoint(mdInitialiseDict_.lookup("startPoint"));

    vector lengthVector = mdInitialiseDict_.lookup("lengthVector");
    vector breadthVector = mdInitialiseDict_.lookup("breadthVector");
    vector normalVector = mdInitialiseDict_.lookup("thicknessVector");

    scalar length = readScalar(mdInitialiseDict_.lookup("length"));
    scalar breadth = readScalar(mdInitialiseDict_.lookup("breadth"));
    scalar thickness = readScalar(mdInitialiseDict_.lookup("thickness"));

    breadthVector /= mag(breadthVector);
    lengthVector /= mag(lengthVector);
    normalVector /= mag(normalVector);

    {
        scalar dotProd = breadthVector & lengthVector;
    
        if(dotProd > SMALL)
        {
            FatalErrorIn("polySimpleLattice::setInitialConfiguration()")
                << "breadthVector: " << breadthVector 
                << " and lengthVector: " << lengthVector
                << " have to be perpendicular to each other. Value of dot product: " 
                << dotProd
                << nl << "in mdInitialiseDict."
                << exit(FatalError);
        }
    }


    {
        scalar dotProd = normalVector & lengthVector;
    
        if(dotProd > SMALL)
        {
            FatalErrorIn("polySimpleLattice::setInitialConfiguration()")
                << "thicknessVector: " << normalVector 
                << " and lengthVector: " << lengthVector
                << " have to be perpendicular to each other. Value of dot product: " 
                << dotProd
                << nl << "in mdInitialiseDict."
                << exit(FatalError);
        }
    }

    if(computeSpacingFromDensity)
    {
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
                    << abort(FatalError);
            }
        }
        else if (mdInitialiseDict_.found("massDensity"))
        {
//             const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));
    
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
                    << abort(FatalError);
            }
        }
        else if (mdInitialiseDict_.found("massDensitySI"))
        {
//             const polyMolecule::constantProperties& cP(molCloud_.constProps(molId));
    
            scalar mass = molCloud_.cP().mass(molId);
    
            Info << "mass: " << mass << endl;
    
            scalar massDensity = readScalar
            (
                mdInitialiseDict_.lookup("massDensitySI")
            );
            
            const reducedUnits& rU =molCloud_.redUnits();

            massDensity /= rU.refMassDensity();
            
            numberDensity = massDensity / mass;
            
            Info << " number density in reduced units = " << numberDensity << endl;
            
            if (massDensity < VSMALL)
            {
                FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                    << "massDensity too small, not filling zone "
                    << abort(FatalError);
            }
        }
        else
        {
            FatalErrorIn("Foam::polyMoleculeCloud::initialiseMolecules")
                << "massDensity or numberDensity not specified " << nl
                << abort(FatalError);
        }
   
        spacing = pow(  (1.0/numberDensity), (1.0/3.0) );
    }

    vector globalPosition = vector::zero;

    label nMolsX = label(length/spacing) + 1;
    label nMolsY = label(breadth/spacing) + 1;
    label nMolsZ = label(thickness/spacing) + 1;

    Info << "Spacing: " << spacing << endl;

    scalar spacingX = spacing;
    scalar spacingY = spacing;
    scalar spacingZ = spacing; 
/*
    if(nMolsX > 1)
    {
        spacingX = length/scalar(nMolsX-1);

        Info<< "nMolsX less than zero: modify spacingX: " << spacingX << endl;
    }

    if(nMolsY > 1)
    {
        spacingY = breadth/scalar(nMolsY-1);

        Info<< "nMolsY less than zero: modify spacingY: " << spacingY << endl;
    }

    if(nMolsZ > 1)
    {
        spacingZ = thickness/scalar(nMolsZ-1);

        Info<< "nMolsZ less than zero: modify spacingZ: " << spacingZ << endl;
    }*/

    if(adjustSpacing)
    {
        spacingX = length/scalar(nMolsX - 1);
        spacingY = breadth/scalar(nMolsY - 1);
        spacingZ = thickness/scalar(nMolsZ - 1);

        Info<< "adjust spacing ---> spacingX: " << spacingX 
            << ", spacingY: " << spacingY
            << ", spacingZ: " << spacingZ
            << endl;
    }





    scalar x = 0.0;
    scalar y = 0.0;
    scalar z = 0.0;

    Info << nl << "building lattice.... " << endl;

    for (label i = 0; i < nMolsX; i++)
    {
        x = i*spacingX;

        for (label j = 0; j < nMolsY; j++)
        {
            y = j*spacingY;

            for (label k = 0; k < nMolsZ; k++)
            {
                z = k*spacingZ;

                globalPosition = startPoint 
                                + lengthVector*x 
                                + breadthVector*y 
                                + normalVector*z;
            
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
                    if(deleteOutOfZone)
                    {
                        if(findIndex(cells, cell) != -1)
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
                    else
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
                else if(!deleteOutOfMesh)
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
//                 else if (deleteOutOfMesh)
//                 {
//                     FatalErrorIn("Foam::polySimpleLattice::insertMolecule()")
//                         << "Position of molecule out of mesh: " << globalPosition
//                         << nl
//                         << abort(FatalError);
//                 }
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


    //delete overlaps, run in serial only
//     bool deleteOverlaps = true;

    if(deleteOverlaps)
    {
        Info << nl << "deleting overlaps " << nl << endl;
//         label initialSize = molCloud_.size();

        DynamicList<polyMolecule*> molsToDel;
        DynamicList<label> molsToDelLabel(0);

        IDLList<polyMolecule>::iterator molI(molCloud_.begin());
        label i = 0;
        label j = 0;
   
        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            IDLList<polyMolecule>::iterator molJ(molCloud_.begin());

            j = 0;

            for
            (
                molJ = molCloud_.begin();
                molJ != molCloud_.end();
                ++molJ
            )
            {

                if(j > i)
                {
                    scalar rIJMag = mag(molI().position() - molJ().position());

                    if(rIJMag < rOv)
                    {
                        if(molJ().id() != molId)
                        {
                            if(findIndex(molsToDelLabel, j) == -1)
                            {
                                polyMolecule* molD = &molJ();
                                molsToDel.append(molD);
                                molsToDelLabel.append(j);
                            }
                        }
                    }
                }

                j++;
            }

            i++;
        }
    
        //molsToDel.shrink();

        label deletedMols = molsToDel.size();

        forAll(molsToDel, m)
        {
            deleteMolecule(*molsToDel[m]);
        }
   
        Info<< tab << " overlapping molecules: " <<  deletedMols
            << endl;
    }
}




} // End namespace Foam

// ************************************************************************* //
