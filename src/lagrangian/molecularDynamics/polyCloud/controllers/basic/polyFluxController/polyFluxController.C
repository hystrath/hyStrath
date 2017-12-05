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

#include "polyFluxController.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyFluxController, 0);

defineRunTimeSelectionTable(polyFluxController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyFluxController::polyFluxController
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
    controllerDict_(dict.subDict("controllerProperties")),
    time_(t),
	regionName_(controllerDict_.lookup("zoneName")),
    regionId_(-1),
    zoneSurfaceArea_(0.0),
    internalFaces_(),
    processorFaces_(),
    control_(true),
    readStateFromFile_(true),
    singleValueController_(false),
    density_(0.0),
    velocity_(vector::zero),
    temperature_(0.0),
    pressure_(0.0),
    strainRate_(tensor::zero),
    tempGradient_(vector::zero),
    fieldController_(false),
    densities_(),
    velocities_(),
    temperatures_(),
    pressures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const faceZoneMesh& faceZones = mesh_.faceZones();
    regionId_ = faceZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyFluxController::polyFluxController()")
        << "Cannot find region (faceZone): " << regionName_ << nl << "in: "
        << t.time().system()/"controllersDict"
        << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    setFacesInfo();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyFluxController> polyFluxController::New
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyFluxControllerName
    (
        dict.lookup("fluxControllerModel")
    );

    Info<< "Selecting fluxController "
         << polyFluxControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyFluxControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyFluxController::New(const dictionary&) : " << endl
            << "    unknown polyFluxController type "
            << polyFluxControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyFluxController>
	(
		cstrIter()(t, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyFluxController::~polyFluxController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void polyFluxController::setFacesInfo()
{
    const labelList& faces = controlZone();

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

        processorFaces_.setSize(processorFaces.size(), -1);

        forAll(processorFaces, f)
        {
            processorFaces_[f] = processorFaces[f];
        }

        label nInternalFaces = faces.size() - processorFaces.size();
        internalFaces_.setSize(nInternalFaces, -1);

        label counter = 0;

        forAll(faces, f)
        {
            const label& faceI = faces[f];

            if(findIndex(processorFaces, faceI) == -1)
            {
                internalFaces_[counter] = faceI;
                counter++;
            }
        }

        forAll(internalFaces_, f)
        {
            const label& faceI = internalFaces_[f];
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }

        // faces on a zone located on a processor cut belong to both processors (hence the 0.5)

        forAll(processorFaces_, f)
        {
            const label& faceI = processorFaces_[f];
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


void polyFluxController::updateFluxControllerProperties
(
    const dictionary& newDict
)
{
    controllerDict_ = newDict.subDict("controllerProperties");

    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can infact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

    if (controllerDict_.found("controlSwitch"))
    {
        control_ = Switch(controllerDict_.lookup("controlSwitch"));
    }

    if (controllerDict_.found("readStateFromFile"))
    {
        readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));
    }
}

const labelList& polyFluxController::controlZone() const
{
    return mesh_.faceZones()[regionId_];
}

label polyFluxController::isFaceOnControlZone(const label& faceI) 
{
    const label f = findIndex(controlZone(), faceI);

    return f;
}

const word& polyFluxController::regionName() const
{
    return regionName_;
}

const scalar& polyFluxController::density() const
{
    return density_;
}

scalar& polyFluxController::density()
{
    return density_;
}

const vector& polyFluxController::velocity() const
{
    return velocity_;
}

vector& polyFluxController::velocity()
{
    return velocity_;
}

const scalar& polyFluxController::temperature() const
{
    return temperature_;
}

scalar& polyFluxController::temperature()
{
    return temperature_;
}

const scalar& polyFluxController::pressure() const
{
    return pressure_;
}

scalar& polyFluxController::pressure()
{
    return pressure_;
}

const tensor& polyFluxController::strainRate() const
{
    return strainRate_;
}

tensor& polyFluxController::strainRate()
{
    return strainRate_;
}

const vector& polyFluxController::tempGradient() const
{
    return tempGradient_;
}

vector& polyFluxController::tempGradient()
{
    return tempGradient_;
}


const scalarField& polyFluxController::densityField() const
{
    return densities_;
}

scalarField& polyFluxController::densityField()
{
    return densities_;
}

const vectorField& polyFluxController::velocityField() const
{
    return velocities_;
}
vectorField& polyFluxController::velocityField()
{
    return velocities_;
}

const scalarField& polyFluxController::temperatureField() const
{
    return temperatures_;
}

scalarField& polyFluxController::temperatureField()
{
    return temperatures_;
}

const scalarField& polyFluxController::pressureField() const
{
    return pressures_;
}

scalarField& polyFluxController::pressureField()
{
    return pressures_;
}


const bool& polyFluxController::singleValueController() const
{
    return singleValueController_;
}

bool& polyFluxController::singleValueController()
{
    return singleValueController_;
}

const bool& polyFluxController::fieldController() const
{
    return fieldController_;
}

bool& polyFluxController::fieldController()
{
    return fieldController_;
}


const bool& polyFluxController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyFluxController::writeInCase() const
{
    return writeInCase_;
}

} // End namespace Foam

// ************************************************************************* //
