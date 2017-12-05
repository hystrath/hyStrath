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

#include "dsmcStateController.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcStateController, 0);

defineRunTimeSelectionTable(dsmcStateController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcStateController::dsmcStateController
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(cloud.mesh()),
    cloud_(cloud),
    rndGen_(cloud_.rndGen()),
    controllerDict_(dict.subDict("controllerProperties")),
    timeDict_(controllerDict_.subDict("timeProperties")),
    time_(t, timeDict_),
    timePeriod_(readScalar(timeDict_.lookup("initialTimePeriod"))), //temp
    initialTime_(time_.time().startTime().value()),
    regionName_(controllerDict_.lookup("zoneName")),
    regionId_(-1),
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
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcStateController::dsmcStateController()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    const scalar& avTimeInterval = time_.averageTimeInterval().deltaT();

    if((timePeriod_ < avTimeInterval) && (timePeriod_ > 0.0))
    {
        timePeriod_ = avTimeInterval;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcStateController> dsmcStateController::New
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcStateControllerName
    (
        dict.lookup("stateControllerModel")
    );

    Info<< "Selecting dsmcStateController "
         << dsmcStateControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcStateControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcStateController::New(const dictionary&) : " << endl
            << "    unknown dsmcStateController type "
            << dsmcStateControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcStateController>
	(
		cstrIter()(t, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcStateController::~dsmcStateController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcStateController::updateTime()
{
    time_++;

    const scalar& t = time_.time().timeOutputValue();
    
    if((t - initialTime_) < timePeriod_)
    {
        time_.controlTimeInterval().endTime() = false;
//         control_ = false;
    }
    else
    {
//         control_ = true;
    }
}

void dsmcStateController::writeTimeData
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
        FatalErrorIn("dsmcStateController::writeTimeData()")
            << "Cannot open file " << timeFile.name()
            << abort(FatalError);
    }
}

void dsmcStateController::writeTimeData
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

void dsmcStateController::writeTimeData
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
        FatalErrorIn("void dsmcStateController::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

void dsmcStateController::writeTimeData
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
        FatalErrorIn("void dsmcStateController::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

void dsmcStateController::updateStateControllerProperties
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

    timeDict_ = controllerDict_.subDict("timeProperties");

    if (timeDict_.found("resetAtOutput"))
    {
        time_.resetFieldsAtOutput() = Switch(timeDict_.lookup("resetAtOutput"));
    }
}

const labelList& dsmcStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

const word& dsmcStateController::regionName() const
{
    return regionName_;
}

const scalar& dsmcStateController::density() const
{
    return density_;
}

scalar& dsmcStateController::density()
{
    return density_;
}

const vector& dsmcStateController::velocity() const
{
    return velocity_;
}

vector& dsmcStateController::velocity()
{
    return velocity_;
}

const scalar& dsmcStateController::temperature() const
{
    return temperature_;
}

scalar& dsmcStateController::temperature()
{
    return temperature_;
}

const scalar& dsmcStateController::pressure() const
{
    return pressure_;
}

scalar& dsmcStateController::pressure()
{
    return pressure_;
}

const tensor& dsmcStateController::strainRate() const
{
    return strainRate_;
}

tensor& dsmcStateController::strainRate()
{
    return strainRate_;
}

const vector& dsmcStateController::tempGradient() const
{
    return tempGradient_;
}

vector& dsmcStateController::tempGradient()
{
    return tempGradient_;
}

const scalarField& dsmcStateController::densityField() const
{
    return densities_;
}

scalarField& dsmcStateController::densityField()
{
    return densities_;
}

const vectorField& dsmcStateController::velocityField() const
{
    return velocities_;
}
vectorField& dsmcStateController::velocityField()
{
    return velocities_;
}

const scalarField& dsmcStateController::temperatureField() const
{
    return temperatures_;
}

scalarField& dsmcStateController::temperatureField()
{
    return temperatures_;
}


const scalarField& dsmcStateController::pressureField() const
{
    return pressures_;
}

scalarField& dsmcStateController::pressureField()
{
    return pressures_;
}


const bool& dsmcStateController::singleValueController() const
{
    return singleValueController_;
}

bool& dsmcStateController::singleValueController()
{
    return singleValueController_;
}

const bool& dsmcStateController::fieldController() const
{
    return fieldController_;
}

bool& dsmcStateController::fieldController()
{
    return fieldController_;
}

const bool& dsmcStateController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcStateController::writeInCase() const
{
    return writeInCase_;
}



scalar dsmcStateController::avReqDensity()
{
    scalar totalDensity = 0.0;

    if(singleValueController_) 
    {
        totalDensity = density_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(densities_, c)
        {
            totalDensity += densities_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalDensity, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalDensity /= scalar(controlCells);
        }
    }

    return totalDensity;
}

vector dsmcStateController::avReqVelocity()
{
    vector totalVel = vector::zero;

    if(singleValueController_) 
    {
        totalVel = velocity_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(velocities_, c)
        {
            totalVel += velocities_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalVel, sumOp<vector>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalVel /= scalar(controlCells);
        }
    }

    return totalVel;
}


scalar dsmcStateController::avReqTemperature()
{
    scalar totalTemp = 0.0;

    if(singleValueController_) 
    {
        totalTemp = temperature_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(temperatures_, c)
        {
            totalTemp += temperatures_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalTemp, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalTemp /= scalar(controlCells);
        }
    }
    return totalTemp;
}

scalar dsmcStateController::avReqPressure()
{
    scalar totalPressure = 0.0;

    if(singleValueController_) 
    {
        totalPressure = pressure_;
    }
    else if(fieldController_)
    {
        label controlCells = controlZone().size();
    
        forAll(pressures_, c)
        {
            totalPressure += pressures_[c];
        }
    
        if (Pstream::parRun())
        {
            reduce(totalPressure, sumOp<scalar>());
            reduce(controlCells, sumOp<label>());
        }
    
        if(controlCells > 0)
        {
            totalPressure /= scalar(controlCells);
        }
    }

    return totalPressure;
}

} // End namespace Foam

// ************************************************************************* //
