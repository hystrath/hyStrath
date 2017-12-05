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

#include "dsmcMassFluxSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMassFluxSurface, 0);

addToRunTimeSelectionTable(dsmcField, dsmcMassFluxSurface, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMassFluxSurface::dsmcMassFluxSurface
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh,cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    regionId_( -1),
    faceZoneName_(propsDict_.lookup("faceZoneName")),
    zoneSurfaceArea_(0.0),
    typeIds_(),
    fluxDirection_(propsDict_.lookup("fluxDirection")),
    molsZone_(0.0),
    massZone_(0.0),
    averagingCounter_(0.0),
    timeIndex_(0),
    molFluxZone_(1, 0.0),
    massFluxZone_(1, 0.0),
    massFlowZone_(1, 0.0),
    averagingAcrossManyRuns_(false)

{

   // choose molecule ids to sample

    const List<word> molecules (propsDict_.lookup("typeIds"));

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcMassFluxSurface::dsmcMassFluxSurface()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = Switch(propsDict_.lookup("averagingAcrossManyRuns"));
        
        // read in stored data from dictionary
        if(averagingAcrossManyRuns_)
        {
            Info << nl << "Averaging across many runs initiated." << nl << endl;

            readIn();
        }         
    }


    // select face zone

    const faceZoneMesh& faceZones = mesh_.faceZones();

    regionId_ = faceZones.findZoneID(faceZoneName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcMassFluxSurface::dsmcMassFluxSurface()")
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
    
    Info << "zoneSurfaceArea_ = " << zoneSurfaceArea_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMassFluxSurface::~dsmcMassFluxSurface()
{}





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcMassFluxSurface::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "massFluxSurface_"+fieldName_+"_"+faceZoneName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("molsZone", molsZone_);
    dict.readIfPresent("massZone", massZone_);    
    
    dict.readIfPresent("averagingCounter", averagingCounter_);
    
//     Info << "Some properties read in: "
//          << "mols = " << mols_[0] 
//          << ", mass = " << mass_[0]
//          << ", averagingCounter = " << averagingCounter_
//          << endl;
}

void dsmcMassFluxSurface::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "massFluxSurface_"+fieldName_+"_"+faceZoneName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("molsZone", molsZone_);
        dict.add("massZone", massZone_);    
    
        dict.add("averagingCounter", averagingCounter_); 
        
        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
        
//         Info<< "Some properties written out: "
//             << "mols = " << mols_[0]
//             << ", mass = " << mass_[0]
//             << ", averagingCounter = " << averagingCounter_
//             << endl;
    }
}

void dsmcMassFluxSurface::createField()
{
}
//- call this function every time-step before the state and flux objects are cleaned
void dsmcMassFluxSurface::calculateField()
{
    averagingCounter_++;

    const List<scalarField>& molIdFlux = cloud_.tracker().parcelIdFlux();
    const List<scalarField>& massIdFlux = cloud_.tracker().massIdFlux();


    scalar molFlux = 0.0;
    scalar massFlux = 0.0;

    const faceZoneMesh& faceZones = mesh_.faceZones();
    const labelList& faces = faceZones[regionId_];

    forAll(faces, f)
    {
        const label& faceI = faces[f];
        vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

        forAll(molIdFlux, id)
        {
            if(findIndex(typeIds_, id) != -1)
            {
                if(cloud_.axisymmetric())
                {
                    const point& fC = cloud_.mesh().faceCentres()[faceI];
                    scalar radius = fC.y();
                    
                    scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                    
                    molFlux += (molIdFlux[id][faceI]*cloud_.nParticle()*RWF*nF) & fluxDirection_;
                    massFlux += (massIdFlux[id][faceI]*cloud_.nParticle()*RWF*nF) & fluxDirection_;
                }
                else
                {
                    molFlux += (molIdFlux[id][faceI]*cloud_.nParticle()*nF) & fluxDirection_;
                    massFlux += (massIdFlux[id][faceI]*cloud_.nParticle()*nF) & fluxDirection_;
                }
            }
        }
    }

    molsZone_ += molFlux;
    massZone_ += massFlux;

    const Time& runTime = time_.time();
    
    // -average measurement and calculate properties
    if(runTime.outputTime())
    {
        scalar molsZone = molsZone_;
        scalar massZone = massZone_;

        if(Pstream::parRun())
        {
            reduce(molsZone, sumOp<scalar>());
            reduce(massZone, sumOp<scalar>());
        }

        scalar averagingTime = averagingCounter_*time_.mdTimeInterval().deltaT();

        molFluxZone_[timeIndex_] = molsZone/(averagingTime*zoneSurfaceArea_);
        massFluxZone_[timeIndex_] = massZone/(averagingTime*zoneSurfaceArea_);
        massFlowZone_[timeIndex_] = massZone/averagingTime;

        if(time_.resetFieldsAtOutput())
        {
            molsZone_ = 0.0;
            massZone_ = 0.0;
            averagingCounter_ = 0.0;
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }

        timeIndex_++;
    }
}



void dsmcMassFluxSurface::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        timeIndex_ = 0;

        if(Pstream::master())
        {
            scalarField timeField(1);
           
            timeField[0] = time_.time().timeOutputValue();

            writeTimeData
            (
                casePath_,
                "faceFlux_"+faceZoneName_+"_"+fieldName_+"_N.xy",
                timeField,
                molFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceFlux_"+faceZoneName_+"_"+fieldName_+"_M.xy",
                timeField,
                massFluxZone_,
                true
            );
            
            writeTimeData
            (
                casePath_,
                "faceMassFlowRate_"+faceZoneName_+"_"+fieldName_+".xy",
                timeField,
                massFlowZone_,
                true
            );
        }
    }
}


// const propertyField& dsmcMassFluxSurface::fields() const
// {
//     return  fields_;
// }

void dsmcMassFluxSurface::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}

} // End namespace Foam

// ************************************************************************* //
