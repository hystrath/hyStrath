/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

Class
    dsmcFluxSurface

\*---------------------------------------------------------------------------*/

#include "dsmcFluxSurface.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcFluxSurface, 0);

addToRunTimeSelectionTable(dsmcField, dsmcFluxSurface, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcFluxSurface::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "fluxSurface_"+fieldName_+"_"+faceZoneName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("particlesAccumulated", particlesAccumulated_);
    dict.readIfPresent("massAccumulated", massAccumulated_);
    dict.readIfPresent("averagingCounter", averagingCounter_);
}

void dsmcFluxSurface::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "fluxSurface_"+fieldName_+"_"+faceZoneName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("particlesAccumulated", particlesAccumulated_);
        dict.add("massAccumulated", massAccumulated_);
        dict.add("averagingCounter", averagingCounter_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

scalar dsmcFluxSurface::calculateFaceZoneArea
(
    const labelList& faces
)
{
    // as outlined in the old dsmcMassFluxSurface implementation special care
    // has to be taken when calculating the faceZone area as both processors
    // will contribute the total area of the face separating them.
    scalar faceZoneArea = 0.0;

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

        const label nInternalFaces = faces.size() - processorFaces.size();
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
            faceZoneArea += mag(mesh_.faceAreas()[faceI]);
        }

        forAll(processorFaces, f)
        {
            const label& faceI = processorFaces[f];
            faceZoneArea += 0.5*mag(mesh_.faceAreas()[faceI]);
        }

        if(Pstream::parRun())
        {
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const label proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << faceZoneArea;
                    }
                }
            }

            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar faceZoneAreaProc;

                    const label proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> faceZoneAreaProc;
                    }
                    faceZoneArea += faceZoneAreaProc;
                }
            }
        }
    }
    else
    {
        forAll(faces, f)
        {
            const label& faceI = faces[f];
            faceZoneArea += mag(mesh_.faceAreas()[faceI]);
        }
    }

    return faceZoneArea;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFluxSurface::dsmcFluxSurface
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    regionId_(-1),
    faceZoneName_(propsDict_.lookup("faceZoneName")),
    faceZoneArea_(0.0),
    typeIds_(),
    fluxDirection_(propsDict_.lookup("fluxDirection")),
    particlesAccumulated_(0.0),
    massAccumulated_(0.0),
    averagingCounter_(0.0),
    timeIndex_(0),
    particleFlux_(1, 0.0),
    particleFlowRate_(1, 0.0),
    massFlux_(1, 0.0),
    massFlowRate_(1, 0.0),
    averagingAcrossManyRuns_(false)
{
   // get list of typeIds that should be sampled
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
            FatalErrorInFunction
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << time_.time().system() << "/fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = Switch(propsDict_.lookup("averagingAcrossManyRuns"));

        // in case averaging across many runs is active, read in infos stored
        // in the previous run.
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
        FatalErrorInFunction
            << "Cannot find region: " << faceZoneName_ << nl << "in: "
            << time_.time().system() << "/fieldPropertiesDict"
            << exit(FatalError);
    }

    // normalize flux direction vector
    fluxDirection_ /= mag(fluxDirection_);

    // calculate total area of all faces in the faceZone
    const labelList& faces = faceZones[regionId_];
    faceZoneArea_ = calculateFaceZoneArea(faces);

    Info << "dsmcFluxSurface total area = " << faceZoneArea_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFluxSurface::~dsmcFluxSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcFluxSurface::createField()
{
}

// call this function every timestep before the state and flux objects are cleaned
void dsmcFluxSurface::calculateField()
{
    averagingCounter_++;

    // note: parcelIdFluxes does already contain RWFs
    const List<scalarField>& parcelIdFluxes = cloud_.tracker().parcelIdFlux();
    const List<scalarField>& massIdFluxes = cloud_.tracker().massIdFlux();
    // TODO: obviously other particle properties like momentum, energy, etc.
    // could also be counted here, cf. dsmcFaceTracker.

    scalar particles = 0.0;
    scalar mass = 0.0;

    const faceZoneMesh& faceZones = mesh_.faceZones();
    const labelList& faces = faceZones[regionId_];

    forAll(faces, f)
    {
        const label& faceI = faces[f];
        const vector nF =
            mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

        const label& cellI = mesh_.faceOwner()[faceI];
        // note: we have to use the dtModel().nParticles(cellI) call here
        // as it must not include the RWFs (as outlined above these are
        // already included in the parcelIdFluxes and massIdFluxes)
        const scalar& nParticles = cloud_.coordSystem().dtModel().nParticles(cellI);

        forAll(parcelIdFluxes, id)
        {
            if(findIndex(typeIds_, id) != -1)
            {
                particles +=
                    (parcelIdFluxes[id][faceI]*nParticles*nF) & fluxDirection_;
                mass +=
                    (massIdFluxes[id][faceI]*nParticles*nF) & fluxDirection_;
            }
        }
    }

    particlesAccumulated_ += particles;
    massAccumulated_ += mass;

    const Time& runTime = time_.time();

    // averaging of the accumulated values at output time
    if(runTime.outputTime())
    {
        scalar particlesAccumulated = particlesAccumulated_;
        scalar massAccumulated = massAccumulated_;

        if(Pstream::parRun())
        {
            reduce(particlesAccumulated, sumOp<scalar>());
            reduce(massAccumulated, sumOp<scalar>());
        }

        const scalar averagingTime =
            averagingCounter_*time_.mdTimeInterval().deltaT();

        // particle based quantities:
        const scalar particleFlowRate = particlesAccumulated / averagingTime;
        particleFlowRate_[timeIndex_] = particleFlowRate;
        particleFlux_[timeIndex_] = particleFlowRate / faceZoneArea_;

        // mass based quantities:
        const scalar massFlowRate = massAccumulated / averagingTime;
        massFlowRate_[timeIndex_] = massFlowRate;
        massFlux_[timeIndex_] = massFlowRate / faceZoneArea_;

        if(time_.resetFieldsAtOutput())
        {
            particlesAccumulated_ = 0.0;
            massAccumulated_ = 0.0;
            averagingCounter_ = 0.0;
        }

        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }

        timeIndex_++;
    }
}

void dsmcFluxSurface::resetField()
{
    // FIXME: This has not been properly tested, therefore warn until this is
    // verified to work correctly.
    WarningInFunction
        << "This functionality has not been tested thoroughly and might not "
        << "behave as intended. Use at your own risk and validate the results."
        << endl;

    // this is called if the mesh has been changed (e.g. dynamic mesh
    // adjustment). Therefore the fields have to be resetted as e.g. the number
    // of cells can change.

    // read stored data
    if(averagingAcrossManyRuns_)
    {
        readIn();
    }

    // check if the faceZone still exists in the new mesh
    const faceZoneMesh& faceZones = mesh_.faceZones();
    regionId_ = faceZones.findZoneID(faceZoneName_);

    if(regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << faceZoneName_ << nl << "in "
            << time_.time().system() << "/fieldPropertiesDict"
            << exit(FatalError);
    }

    // recalculate total area of all faces in the faceZone, but check that the
    // total area did not change after the mesh alteration. Therefore save it
    // before recalculation.
    const scalar faceZoneAreaBeforeReset = faceZoneArea_;

    // recalculate the total area of all faces in the faceZone from scratch
    const labelList& faces = faceZones[regionId_];
    faceZoneArea_ = calculateFaceZoneArea(faces);

    // check that faceZoneArea did not change (significantly)
    if (notEqual(faceZoneArea_, faceZoneAreaBeforeReset))
    {
        FatalErrorInFunction
        << "dsmcFluxSurface total area after reset (" << faceZoneArea_ << ") "
        << "is not equal to total area before reset ("
        << faceZoneAreaBeforeReset << ")!"
        << exit(FatalError);
    }

    Info << "dsmcFluxSurface total area = " << faceZoneArea_ << nl << endl;
}

void dsmcFluxSurface::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        timeIndex_ = 0;

        if(Pstream::master())
        {
            scalarField timeField(1);

            timeField[0] = time_.time().timeOutputValue();

            // .xy extension chosen to be in line with typical OpenFOAM
            // sampling
            writeTimeData
            (
                casePath_,
                "faceParticleFlux_" + faceZoneName_ + "_" + fieldName_ + ".xy",
                timeField,
                particleFlux_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceParticleFlowRate_" + faceZoneName_ + "_" + fieldName_
                    + ".xy",
                timeField,
                particleFlowRate_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceMassFlux_" + faceZoneName_ + "_" + fieldName_ + ".xy",
                timeField,
                massFlux_,
                true
            );

            writeTimeData
            (
                casePath_,
                "faceMassFlowRate_" + faceZoneName_ + "_" + fieldName_ + ".xy",
                timeField,
                massFlowRate_,
                true
            );
        }
    }
}

void dsmcFluxSurface::updateProperties(const dictionary& newDict)
{
    updateBasicFieldProperties(newDict);
}

// Note: This is pretty hacky. These functions are called by the fields_ object
// in dsmcCloud. The fields_ object is an object of class dsmcFieldProperties.
// This object is responsible for creating, updating, etc. of all the fields in
// fieldPropertiesDict. If fields_.translationalT(cellI) is called it returns
// the translationalT of the FIRST field specified in fieldPropertiesDict.
// Therefore this field should be of type dsmcVolFields. Hence these functions
// will crash and inform the user in case this ordering is not respected.
scalar dsmcFluxSurface::translationalT(const label cellI)
{
    FatalErrorInFunction
        << "Please make sure that the first field in " << time_.time().system()
        << "/fieldPropertiesDict is of type dsmcVolFields."
        << exit(FatalError);

    return -1;
}

scalar dsmcFluxSurface::overallT(const label cellI)
{
    FatalErrorInFunction
        << "Please make sure that the first field in " << time_.time().system()
        << "/fieldPropertiesDict is of type dsmcVolFields."
        << exit(FatalError);

    return -1;
}

} // End namespace Foam

// ************************************************************************* //
