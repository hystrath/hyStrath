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

#include "dsmcFluxController.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcFluxController, 0);

defineRunTimeSelectionTable(dsmcFluxController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFluxController::dsmcFluxController
(
    Time& t,
    dsmcCloud& cloud,
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
        FatalErrorIn("dsmcFluxController::dsmcFluxController()")
        << "Cannot find region (faceZone): " << regionName_ << nl << "in: "
        << t.time().system()/"controllersDict"
        << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    setFacesInfo();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcFluxController> dsmcFluxController::New
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcFluxControllerName
    (
        dict.lookup("fluxControllerModel")
    );

    Info<< "Selecting fluxController "
         << dsmcFluxControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcFluxControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcFluxController::New(const dictionary&) : " << endl
            << "    unknown dsmcFluxController type "
            << dsmcFluxControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcFluxController>
	(
		cstrIter()(t, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFluxController::~dsmcFluxController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void dsmcFluxController::updateTime()
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


void dsmcFluxController::setFacesInfo()
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

//     Pout << "ERROR 1" << endl;
}

void dsmcFluxController::writeTimeData
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

void dsmcFluxController::writeTimeData
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

void dsmcFluxController::writeTimeData
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

void dsmcFluxController::writeTimeData
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

void dsmcFluxController::updateTime()
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

void dsmcFluxController::updateFluxControllerProperties
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

const labelList& dsmcFluxController::controlZone() const
{
    return mesh_.faceZones()[regionId_];
}

label dsmcFluxController::isFaceOnControlZone(const label& faceI) 
{
    const label f = findIndex(controlZone(), faceI);

    return f;
}

const word& dsmcFluxController::regionName() const
{
    return regionName_;
}

const scalar& dsmcFluxController::density() const
{
    return density_;
}

scalar& dsmcFluxController::density()
{
    return density_;
}

const vector& dsmcFluxController::velocity() const
{
    return velocity_;
}

vector& dsmcFluxController::velocity()
{
    return velocity_;
}

const scalar& dsmcFluxController::temperature() const
{
    return temperature_;
}

scalar& dsmcFluxController::temperature()
{
    return temperature_;
}

const scalar& dsmcFluxController::pressure() const
{
    return pressure_;
}

scalar& dsmcFluxController::pressure()
{
    return pressure_;
}

const tensor& dsmcFluxController::strainRate() const
{
    return strainRate_;
}

tensor& dsmcFluxController::strainRate()
{
    return strainRate_;
}

const vector& dsmcFluxController::tempGradient() const
{
    return tempGradient_;
}

vector& dsmcFluxController::tempGradient()
{
    return tempGradient_;
}


const scalarField& dsmcFluxController::densityField() const
{
    return densities_;
}

scalarField& dsmcFluxController::densityField()
{
    return densities_;
}

const vectorField& dsmcFluxController::velocityField() const
{
    return velocities_;
}
vectorField& dsmcFluxController::velocityField()
{
    return velocities_;
}

const scalarField& dsmcFluxController::temperatureField() const
{
    return temperatures_;
}

scalarField& dsmcFluxController::temperatureField()
{
    return temperatures_;
}

const scalarField& dsmcFluxController::pressureField() const
{
    return pressures_;
}

scalarField& dsmcFluxController::pressureField()
{
    return pressures_;
}


const bool& dsmcFluxController::singleValueController() const
{
    return singleValueController_;
}

bool& dsmcFluxController::singleValueController()
{
    return singleValueController_;
}

const bool& dsmcFluxController::fieldController() const
{
    return fieldController_;
}

bool& dsmcFluxController::fieldController()
{
    return fieldController_;
}


const bool& dsmcFluxController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcFluxController::writeInCase() const
{
    return writeInCase_;
}


// const scalar dsmcFluxController::avReqDensity() const
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
// const vector dsmcFluxController::avReqVelocity() const
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
// const scalar dsmcFluxController::avReqTemperature() const
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
