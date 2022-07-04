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

Description

\*---------------------------------------------------------------------------*/

#include "pdFluxController.H"
#include "IFstream.H"
#include "graph.H"
#include "pdCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdFluxController, 0);

defineRunTimeSelectionTable(pdFluxController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdFluxController::pdFluxController
(
    Time& t,
    pdCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    rndGen_(clock::getTime()),
    controllerDict_(dict.subDict("controllerProperties")),
    timeDict_(controllerDict_.subDict("timeProperties")),
    time_(t, timeDict_),
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
        FatalErrorIn("pdFluxController::pdFluxController()")
        << "Cannot find region (faceZone): " << regionName_ << nl << "in: "
        << t.time().system()/"controllersDict"
        << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    setFacesInfo();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pdFluxController> pdFluxController::New
(
    Time& t,
    pdCloud& cloud,
    const dictionary& dict
)
{
    word pdFluxControllerName
    (
        dict.lookup("fluxControllerModel")
    );

    Info<< "Selecting fluxController "
         << pdFluxControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdFluxControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdFluxController::New(const dictionary&) : " << endl
            << "    unknown pdFluxController type "
            << pdFluxControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pdFluxController>
	(
		cstrIter()(t, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdFluxController::~pdFluxController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void pdFluxController::updateTime()
// {
//     time_++;
//
//     const scalar& t = time_.time().timeOutputValue();
//
//     if((t - initialTime_) < timePeriod_)
//     {
//         time_.controlTimeInterval().endTime() = false;
// //         control_ = false;
//     }
//     else
//     {
// //         control_ = true;
//     }
// }


void pdFluxController::setFacesInfo()
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

        processorFaces.shrink();

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

//         Pout << "processorFaces: " << processorFaces_ << endl;
//         Pout << "internalFaces: " << internalFaces_ << endl;

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
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const label proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << zoneSurfaceArea_;
                    }
                }
            }

            //- receiving
            for (label p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalar zoneSurfaceAreaProc;

                    const label proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
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

//     Pout << "ERROR 1" << endl;
}

void pdFluxController::writeTimeData
(
    const fileName& pathName,
    const word& fileName,
    const List< Pair<scalar> >& data
)
{
    OFstream timeFile(pathName/fileName+".raw");

    if (timeFile.good())
    {
        timeFile << data << endl;
    }
    else
    {
        FatalErrorIn("stateController::writeTimeData()")
            << "Cannot open file " << timeFile.name()
            << abort(FatalError);
    }
}

void pdFluxController::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData
)
{
    fileName writeFile(pathName/nameFile);

    graph outputGraph("title", "x", "y", xData, yData);

    outputGraph.write(writeFile, "raw");
}

void pdFluxController::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if(file.good())
    {
        forAll(yData, n)
        {
            file
                << xData[n] << "\t"
                << yData[n].x() << "\t" << yData[n].y()
                << "\t" << yData[n].z()
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void stateController::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

void pdFluxController::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const tensorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if(file.good())
    {
        forAll(yData, n)
        {
            file
                << xData[n] << "\t"
                << yData[n].xx() << "\t" << yData[n].xy() << "\t" << yData[n].xz() << "\t"
                << yData[n].yx() << "\t" << yData[n].yy() << "\t" << yData[n].yz() << "\t"
                << yData[n].zx() << "\t" << yData[n].zy() << "\t" << yData[n].zz()
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void stateController::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

void pdFluxController::updateTime()
{
    time_++;

//     const scalar& t = time_.time().timeOutputValue();
//
//     if((t - initialTime_) < timePeriod_)
//     {
//         time_.controlTimeInterval().endTime() = false;
// //         control_ = false;
//     }
//     else
//     {
// //         control_ = true;
//     }
}

void pdFluxController::updateFluxControllerProperties
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

const labelList& pdFluxController::controlZone() const
{
    return mesh_.faceZones()[regionId_];
}

label pdFluxController::isFaceOnControlZone(const label& faceI)
{
    const label f = findIndex(controlZone(), faceI);

    return f;
}

const word& pdFluxController::regionName() const
{
    return regionName_;
}

const scalar& pdFluxController::density() const
{
    return density_;
}

scalar& pdFluxController::density()
{
    return density_;
}

const vector& pdFluxController::velocity() const
{
    return velocity_;
}

vector& pdFluxController::velocity()
{
    return velocity_;
}

const scalar& pdFluxController::temperature() const
{
    return temperature_;
}

scalar& pdFluxController::temperature()
{
    return temperature_;
}

const scalar& pdFluxController::pressure() const
{
    return pressure_;
}

scalar& pdFluxController::pressure()
{
    return pressure_;
}

const tensor& pdFluxController::strainRate() const
{
    return strainRate_;
}

tensor& pdFluxController::strainRate()
{
    return strainRate_;
}

const vector& pdFluxController::tempGradient() const
{
    return tempGradient_;
}

vector& pdFluxController::tempGradient()
{
    return tempGradient_;
}


const scalarField& pdFluxController::densityField() const
{
    return densities_;
}

scalarField& pdFluxController::densityField()
{
    return densities_;
}

const vectorField& pdFluxController::velocityField() const
{
    return velocities_;
}
vectorField& pdFluxController::velocityField()
{
    return velocities_;
}

const scalarField& pdFluxController::temperatureField() const
{
    return temperatures_;
}

scalarField& pdFluxController::temperatureField()
{
    return temperatures_;
}

const scalarField& pdFluxController::pressureField() const
{
    return pressures_;
}

scalarField& pdFluxController::pressureField()
{
    return pressures_;
}


const bool& pdFluxController::singleValueController() const
{
    return singleValueController_;
}

bool& pdFluxController::singleValueController()
{
    return singleValueController_;
}

const bool& pdFluxController::fieldController() const
{
    return fieldController_;
}

bool& pdFluxController::fieldController()
{
    return fieldController_;
}


const bool& pdFluxController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& pdFluxController::writeInCase() const
{
    return writeInCase_;
}


// const scalar pdFluxController::avReqDensity() const
// {
//     scalar totalDensity = 0.0;
//
//     forAll(densities_, c)
//     {
//         totalDensity += densities_[c];
//     }
//
//     if(cells_.size() > 0)
//     {
//         totalDensity /= scalar(cells_.size());
//     }
//
//     return totalDensity;
// }
//
// const vector pdFluxController::avReqVelocity() const
// {
//     vector totalVel = vector::zero;
//
//     forAll(velocities_, c)
//     {
//         totalVel += velocities_[c];
//     }
//
//     if(cells_.size() > 0)
//     {
//         totalVel /= scalar(cells_.size());
//     }
//
//     return totalVel;
// }
//
// const scalar pdFluxController::avReqTemperature() const
// {
//     scalar totalTemp = 0.0;
//
//     forAll(densities_, c)
//     {
//         totalTemp += temperatures_[c];
//     }
//
//     if(cells_.size() > 0)
//     {
//         totalTemp /= scalar(cells_.size());
//     }
//
//     return totalTemp;
// }



} // End namespace Foam

// ************************************************************************* //
