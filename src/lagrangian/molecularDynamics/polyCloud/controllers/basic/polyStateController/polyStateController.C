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

#include "polyStateController.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"


namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyStateController, 0);

defineRunTimeSelectionTable(polyStateController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyStateController::polyStateController
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(molCloud.mesh())),
    molCloud_(molCloud),
//     controllerDict_(dict.subDict("controllerProperties")),
	time_(t), 
    regionName_(dict.lookup("zoneName")),
    regionId_(-1),
    control_(true),
    controlInterForces_(false),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyStateController::polyStateController()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }
    
    if(dict.found("control"))
    {    
        control_ = Switch(dict.lookup("control"));
    }
//     readStateFromFile_ = Switch(dict.lookup("readStateFromFile"));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyStateController> polyStateController::New
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyStateControllerName
    (
        dict.lookup("stateControllerModel")
    );

    Info<< "Selecting polyStateController "
         << polyStateControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyStateControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyStateController::New(const dictionary&) : " << endl
            << "    unknown polyStateController type "
            << polyStateControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyStateController>
	(
		cstrIter()(t, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyStateController::~polyStateController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void polyStateController::updateStateControllerProperties
(
    const dictionary& newDict
)
{
//     controllerDict_ = newDict.subDict("controllerProperties");

    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can infact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

    if (newDict.found("control"))
    {
        control_ = Switch(newDict.lookup("control"));
    }
/*
    if (newDict.found("readStateFromFile"))
    {
        readStateFromFile_ = Switch(newDict.lookup("readStateFromFile"));
    }*/
}

const labelList& polyStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

const word& polyStateController::regionName() const
{
    return regionName_;
}

const bool& polyStateController::controlInterForces() const
{
    return controlInterForces_;
}

bool& polyStateController::controlInterForces()
{
    return controlInterForces_;
}


const bool& polyStateController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyStateController::writeInCase() const
{
    return writeInCase_;
}



} // End namespace Foam

// ************************************************************************* //
