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

#include "polyMassFluxSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMassFluxSurface, 0);
addToRunTimeSelectionTable(polyField, polyMassFluxSurface, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMassFluxSurface::polyMassFluxSurface
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
    molsZone_(0.0),
    massZone_(0.0),
    absMomZone_(0.0),
    momZone_(vector::zero),
    averagingTime_(0.0),
    writeInterval_(readScalar(t.controlDict().lookup("writeInterval"))),
    resetAtOutput_(false),
    molFluxZone_(1, 0.0),
    massFluxZone_(1, 0.0),
    absMomFluxZone_(1, 0.0),
    momFluxZone_(1, 0.0),
    computeErrorBars_(false),
    nSteps_(0),
    stepCounter_(0),
    averagingCounter_(0.0),
    meanMassFlux_(0.0),
    stdTerm_(0.0)
{
    averagingAcrossManyRuns_ = false;

    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = Switch(propsDict_.lookup("averagingAcrossManyRuns"));
    }

    if (propsDict_.found("computeErrorBars"))
    {
        computeErrorBars_ = Switch(propsDict_.lookup("computeErrorBars"));

        if(computeErrorBars_)
        {
            Info << " Mass-flux: computing error bars" << nl << endl;

            nSteps_ = readLabel(propsDict_.lookup("nSteps"));
            meanMassFlux_= readScalar(propsDict_.lookup("meanMassFlux"));
        }
    }
    
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    // read in stored data from dictionary
    if(averagingAcrossManyRuns_)
    {
        Info << " Averaging across many runs. Reading from dictionary:" << endl;

        pathName_ = time_.time().path()/"storage";
        nameFile_ = "massFluxData";

        if( !isDir(pathName_) )
        {
            mkDir(pathName_);

            Info << nl << "Storage not found!"  << nl << endl;
            Info << ".... creating"  << nl << endl;
        }

        bool fileFound = readFromStorage();

        if(!fileFound)
        {
            Info << nl << "File not found: " << nameFile_ << nl << endl;
            Info << ".... creating"  << nl << endl;
            writeToStorage();

            Info << "setting properties to default values. " << endl;
        }
        else
        {
            Pout<< "Properties read-in are: mols = " << molsZone_ << ", mass = " << massZone_
                << ", averagingTime = " << averagingTime_
                << endl;
        }
       
    }

    if (propsDict_.found("resetAtOutput"))
    {
        resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));
    }

    // select face zone

    const faceZoneMesh& faceZones = mesh_.faceZones();

    regionId_ = faceZones.findZoneID(faceZoneName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyMassFluxSurface::polyMassFluxSurface()")
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
        
        //processorFaces.shrink();

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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMassFluxSurface::~polyMassFluxSurface()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMassFluxSurface::createField()
{}

//- call this function every time-step before the state and flux objects are cleaned
void polyMassFluxSurface::calculateField()
{
    const List<scalarField>& molIdFlux = molCloud_.tracker().molIdFlux();
    const List<scalarField>& massIdFlux = molCloud_.tracker().massIdFlux();
    const List<vectorField>& momIdFlux = molCloud_.tracker().momIdFlux();
    const List<scalarField>& absMomIdFlux = molCloud_.tracker().absMomIdFlux();

    scalar molFlux = 0.0;
    scalar massFlux = 0.0;
    scalar absMomFlux = 0.0;
    vector momFlux = vector::zero;

    const faceZoneMesh& faceZones = mesh_.faceZones();
    const labelList& faces = faceZones[regionId_];

    forAll(faces, f)
    {
        const label& faceI = faces[f];
        vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

        forAll(molIdFlux, id)
        {
            if(findIndex(molIds_, id) != -1)
            {
                molFlux += (molIdFlux[id][faceI]*nF) & fluxDirection_;
                massFlux += (massIdFlux[id][faceI]*nF) & fluxDirection_;
                momFlux += momIdFlux[id][faceI];
                absMomFlux += absMomIdFlux[id][faceI];
            }
        }
    }

    molsZone_ += molFlux;
    massZone_ += massFlux;
    momZone_ += momFlux;
    absMomZone_ += absMomFlux;

    if(computeErrorBars_)
    {
        stepCounter_++;

        if(stepCounter_ >= nSteps_)
        {
            stepCounter_ = 0;

            averagingCounter_+=1.0;

            if(Pstream::parRun())
            {
                reduce(massFlux, sumOp<scalar>());
            }
    
            scalar flux = massFlux/time_.deltaT().value(); 
    
            stdTerm_ += (flux - meanMassFlux_)*(flux - meanMassFlux_);
        }
    }


    const Time& runTime = time_;

    // -average measurement and calculate properties
    if(runTime.outputTime())
    {
        scalar molsZone = molsZone_;
        scalar massZone = massZone_;
        scalar absMomZone = absMomZone_;
        vector momZone = momZone_;

        if(Pstream::parRun())
        {
            reduce(molsZone, sumOp<scalar>());
            reduce(massZone, sumOp<scalar>());
            reduce(absMomZone, sumOp<scalar>());
            reduce(momZone, sumOp<vector>());
        }

        averagingTime_ += writeInterval_;

        molFluxZone_[0] = molsZone/averagingTime_;
        massFluxZone_[0] = massZone/averagingTime_;
        absMomFluxZone_[0] = absMomZone/averagingTime_;
        momFluxZone_[0] = (momZone & fluxDirection_)/averagingTime_;

        if(resetAtOutput_)
        {
            averagingTime_ = 0.0;
            molsZone_ = 0.0;
            massZone_ = 0.0;
            absMomZone_ = 0.0;
            momZone_ = vector::zero;
        }

        if(averagingAcrossManyRuns_)
        {
            writeToStorage();
        }
    }
}

void polyMassFluxSurface::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
            file << molsZone_ << endl;
            file << massZone_ << endl;
            file << averagingTime_ << endl;
    }
    else
    {
        FatalErrorIn("void polyMassFluxSurface::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool polyMassFluxSurface::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar mols(readScalar(file));
        scalar mass(readScalar(file));
        scalar averagingTime(readScalar(file));

        molsZone_ = mols;
        massZone_ = mass;
        averagingTime_ = averagingTime;
    }

    return goodFile;
    
}

void polyMassFluxSurface::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());

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
                "faceFlux_"+faceZoneName_+"_"+fieldName_+"_abs_mom.xy",
                timeField,
                absMomFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceFlux_"+faceZoneName_+"_"+fieldName_+"_mom.xy",
                timeField,
                momFluxZone_,
                true
            );

            if(computeErrorBars_)
            {
                scalar sigma = sqrt((stdTerm_)/scalar(averagingCounter_));
    
                scalarField errorBarField(1, 0.0);
    
                errorBarField[0] = sigma/sqrt(scalar(averagingCounter_));
    
                writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneName_+"_"+fieldName_+"_errorBars.xy",
                    timeField,
                    errorBarField,
                    true
                );
            }

            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneName_+"_"+fieldName_+"_N_SI.xy",
                    timeField*rU.refTime(),
                    molFluxZone_*rU.refMolFlux(),
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "faceFlux_"+faceZoneName_+"_"+fieldName_+"_M_SI.xy",
                    timeField*rU.refTime(),
                    massFluxZone_*rU.refMassFlux(),
                    true
                );
            }
        }
    }
}

void polyMassFluxSurface::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyMassFluxSurface::measureDuringForceComputationSite
(   
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyMassFluxSurface::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
