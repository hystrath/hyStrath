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

#include "polyMassFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMassFlux, 0);

addToRunTimeSelectionTable(polyField, polyMassFlux, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMassFlux::polyMassFlux
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),

    fieldName_(propsDict_.lookup("fieldName")),
    fields_(t, mesh, "dummy"),
    regionId_( -1),
    faceZoneName_(propsDict_.lookup("faceZoneName")),
    zoneSurfaceArea_(0.0),
    molIds_(),
    fluxDirection_(propsDict_.lookup("fluxDirection")),
    cumulativeFlux_(0.0)
{
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    // select face zone

    const faceZoneMesh& faceZones = mesh_.faceZones();

    regionId_ = faceZones.findZoneID(faceZoneName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyMassFlux::polyMassFlux()")
            << "Cannot find region: " << faceZoneName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }



    fluxDirection_ /= mag(fluxDirection_);

    // find total surface area
    const labelList& faces = faceZones[regionId_];

    if(Pstream::parRun())
    {
        DynamicList<label> processorFaces(0);

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                for (label p = 0; p < patch.size(); p++)
                {
                    label patchFaceI = p + patch.start();
                    label faceId = findIndex (faces, patchFaceI);

                    if(faceId != -1)
                    {
                        processorFaces.append(patchFaceI);
                    }
                }
            }
        }
        
        processorFaces.shrink();

        label nInternalFaces = faces.size() - processorFaces.size();
           
        List<label> internalFaces(nInternalFaces, 0);

        label counter = 0;

        forAll(faces, f)
        {
            const label& faceI = faces[f];

            if(findIndex(processorFaces, faceI) == -1)
            {
                internalFaces[counter] = faceI;
                counter++;
            }
        }

//         Pout << "processorFaces: " << processorFaces_ << endl;
//         Pout << "internalFaces: " << internalFaces_ << endl;

        forAll(internalFaces, f)
        {
            const label& faceI = internalFaces[f];
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }


    
        forAll(processorFaces, f)
        {
            const label& faceI = processorFaces[f];
            zoneSurfaceArea_ += 0.5*mag(mesh_.faceAreas()[faceI]);           
        }

        if(Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << zoneSurfaceArea_;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar zoneSurfaceAreaProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> zoneSurfaceAreaProc;
                    }
        
                    zoneSurfaceArea_ += zoneSurfaceAreaProc;
                }
            }
        }
    }
    else
    {
        forAll(faces, f)
        {
            const label& faceI = faces[f];

            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }
    }

//     Info << "surface area " << zoneSurfaceArea_ << endl;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMassFlux::~polyMassFlux()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMassFlux::createField()
{
}
//- call this function every time-step before the state and flux objects are cleaned
void polyMassFlux::calculateField()
{
    const List<scalarField>& massIdFlux = molCloud_.tracker().massIdFlux();

    scalar massFlux = 0.0;


    const faceZoneMesh& faceZones = mesh_.faceZones();
    const labelList& faces = faceZones[regionId_];

    forAll(faces, f)
    {
        const label& faceI = faces[f];
        vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

        forAll(massIdFlux, id)
        {
            if(findIndex(molIds_, id) != -1)
            {
                massFlux += (massIdFlux[id][faceI]*nF) & fluxDirection_;
            }
        }
    }

    if(Pstream::parRun())
    {
        reduce(massFlux, sumOp<scalar>());
    }
    
    const scalar& deltaT = time_.time().deltaT().value();
    massFlux /= deltaT;  
    cumulativeFlux_ += massFlux;
    massFlux_.append(massFlux);

    Info << fieldName_ << " - mass flow rate = " << cumulativeFlux_ << endl;
}
    
void polyMassFlux::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            massFlux_.shrink();
            scalarField timeField (massFlux_.size(), 0.0);
            scalarField massFlux (massFlux_.size(), 0.0);
            
            massFlux.transfer(massFlux_);
            massFlux_.clear();

            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeField, i)
            {
                timeField[timeField.size()-i-1]=runTime.timeOutputValue()-(deltaT*i);
            }
            
            writeTimeData
            (
                casePath_,
                "faceFlux_"+faceZoneName_+"_"+fieldName_+"_M.xy",
                timeField,
                massFlux,
                true
            );
                 
            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneName_+"_"+fieldName_+"_M_SI.xy",
                    timeField*rU.refTime(),
                    massFlux*rU.refMassFlux(),
                    true
                );
            }
        }
    }
}

void polyMassFlux::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyMassFlux::measureDuringForceComputationSite
(   
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}

const propertyField& polyMassFlux::fields() const
{
    return  fields_;
}



} // End namespace Foam

// ************************************************************************* //
