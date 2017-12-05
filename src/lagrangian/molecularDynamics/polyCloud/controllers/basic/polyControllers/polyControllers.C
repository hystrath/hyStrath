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

#include "polyControllers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor (for all other md constructors)
polyControllers::polyControllers
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    mesh_(mesh),
    polyControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    nFluxControllers_(0),
    stateControllersList_(),
    sCNames_(),
    sCIds_(),
    sCFixedPathNames_(),
    stateControllers_(),
    fluxControllersList_(),
    fCNames_(),
    fCIds_(),
    fCFixedPathNames_(),
    fluxControllers_(),
    controllersDuringForceComp_()
{}

//- Constructor for mdFOAM
polyControllers::polyControllers
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    time_(t),
    mesh_(mesh),
    polyControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    nFluxControllers_(0),
	stateControllersList_(polyControllersDict_.lookup("polyStateControllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
	stateControllers_(stateControllersList_.size()),
    fluxControllersList_(polyControllersDict_.lookup("polyFluxControllers")),
    fCNames_(fluxControllersList_.size()),
    fCIds_(fluxControllersList_.size()),
    fCFixedPathNames_(fluxControllersList_.size()),
	fluxControllers_(fluxControllersList_.size()),
    controllersDuringForceComp_()
{

    Info << nl << "Creating polyControllers" << nl << endl;

    //- state polyControllers

    DynamicList<label> controllersDuringForceComp(0);

    if(stateControllers_.size() > 0 )
    {
        forAll(stateControllers_, sC)
        {
            const entry& polyControllersI = stateControllersList_[sC];
            const dictionary& polyControllersIDict = polyControllersI.dict();
    
            stateControllers_[sC] = autoPtr<polyStateController>
            (
                polyStateController::New(time_, molCloud, polyControllersIDict)
            );
    
            sCNames_[sC] = stateControllers_[sC]->type();
            sCIds_[sC] = sC;
    
            nStateControllers_++;

            if(stateControllers_[sC]->controlInterForces())
            {
                controllersDuringForceComp.append(sC);
            }
        }
    }

    //controllersDuringForceComp_.transfer(controllersDuringForceComp.shrink());
    controllersDuringForceComp_.transfer(controllersDuringForceComp);

    //- flux polyControllers

    if(fluxControllers_.size() > 0 )
    {
        forAll(fluxControllers_, fC)
        {
            const entry& polyControllersI = fluxControllersList_[fC];
    
            const dictionary& polyControllersIDict = polyControllersI.dict();
    
            fluxControllers_[fC] = autoPtr<polyFluxController>
            (
                polyFluxController::New(time_, molCloud, polyControllersIDict)
            );
    
            fCNames_[fC] = fluxControllers_[fC]->type();
            fCIds_[fC] = fC;
    
            nFluxControllers_++;
        }
    }

    // creating directories for state controllers
    if(nStateControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }

        // directory: case/controllers/poly
        fileName polyControllersPath(controllersPath/"poly");

        if(isDir(polyControllersPath) )
        {
            rmDir(polyControllersPath);
        }

        mkDir(polyControllersPath);

        // directory: case/controllers/poly/stateControllers
        fileName stateControllersPath(polyControllersPath/"stateControllers");
    
        if (!isDir(stateControllersPath))
        {
            mkDir(stateControllersPath);    
        }

        forAll(stateControllers_, sC)
        {
            if(stateControllers_[sC]->writeInCase())
            {
                // directory: case/controllers/poly/stateControllers/<stateControllerModel>
                fileName stateControllerPath(stateControllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);    
                }
    
                const word& regionName = stateControllers_[sC]->regionName();

                // directory: case/controllers/poly/stateControllers/<stateControllerModel>/<cellZoneName>    
                fileName zonePath(stateControllerPath/regionName);
   
                if (!isDir(zonePath))
                {
                    mkDir(zonePath);    
                }
    
                sCFixedPathNames_[sC] = zonePath;
            }
        }
    }

    // creating directories for flux controllers
    if(nFluxControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }

        // directory: case/controllers/poly
        fileName polyControllersPath(time_.path()/"poly");

        if( !isDir(polyControllersPath) )
        {
            mkDir(polyControllersPath);
        }

        // directory: case/controllers/poly/fluxControllers
        fileName fluxControllersPath(polyControllersPath/"fluxControllers");
    
        if (!isDir(fluxControllersPath))
        {
            mkDir(fluxControllersPath);    
        }

        forAll(fluxControllers_, fC)
        {
            if(fluxControllers_[fC]->writeInCase())
            {
                // directory: case/controllers/poly/fluxControllers/<fluxControllerModel>
                fileName fluxControllerPath(fluxControllersPath/fCNames_[fC]);
    
                if (!isDir(fluxControllerPath))
                {
                    mkDir(fluxControllerPath);    
                }

                const word& regionName = fluxControllers_[fC]->regionName();
    
                // directory: case/controllers/poly/fluxControllers/<fluxControllerModel>/<faceZoneName>
                fileName zonePath(fluxControllerPath/regionName);
    
                if (!isDir(zonePath))
                {
                    mkDir(zonePath);    
                }

                fCFixedPathNames_[fC] = zonePath;
            }
        }
    }
}

polyControllers::~polyControllers()
{}

//- initial configuration
//- call this function after the polyMoleculeCloud is completely initialised
void polyControllers::initialConfig()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->initialConfiguration();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->initialConfiguration();
    }
}

void polyControllers::controlVelocitiesI()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeVelocityI();
    }
}

void polyControllers::controlBeforeMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeMove();
    }    
}


void polyControllers::controlBeforeForces()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeForces();
    }
}


void polyControllers::controlDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{
    forAll(controllersDuringForceComp_, n)
    {
        const label& sC = controllersDuringForceComp_[n];
        stateControllers_[sC]->controlDuringForces(molI, molJ);
    }
}

//- control molecular state -- call this after the intermolecular force calulation
void polyControllers::controlAfterForces()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlAfterForces();
    }
}

void polyControllers::controlVelocitiesII()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlAfterVelocityII();
    }
}

//- calculate properties -- call this at the end of the MD time-step.
void polyControllers::calculateStateProps()
{
    forAll(stateControllers_, sC)
    {
//         Info << "error: " << sCNames_[sC] << endl;
        stateControllers_[sC]->calculateProperties();
    }

    forAll(fluxControllers_, fC)
    {
//         Info << "error: " << sCNames_[sC] << endl;
        fluxControllers_[fC]->calculateProperties();
    }
}


//- output -- call this function at the end of the MD time-step
void polyControllers::outputStateResults() 
{
    const Time& runTime = time_;

    if(runTime.outputTime())
    {
        // -- creating a set of directories in the current time directory
        {
            List<fileName> timePathNames(sCFixedPathNames_.size()); 
    
            if(nStateControllers_ > 0)
            {
                if(Pstream::master())
                {
                    // directory: case/<timeDir>/uniform
                    fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
                
                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }
        
        
                    if(stateControllers_.size() > 0)
                    {
                        // directory: case/<timeDir>/uniform/controllers
                        fileName controllersTimePath(uniformTimePath/"controllers");
    
                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }
    
                        // directory: case/<timeDir>/uniform/controllers/poly
                        fileName polyTimePath(controllersTimePath/"poly");
                    
                        if (!isDir(polyTimePath))
                        {
                            mkDir(polyTimePath);    
                        }
    
                        // directory: case/<timeDir>/uniform/controllers/poly/
                        fileName polyStateControllersTimePath(polyTimePath/"stateControllers");
                    
                        if (!isDir(polyStateControllersTimePath))
                        {
                            mkDir(polyStateControllersTimePath);    
                        }
    
                        forAll(stateControllers_, sC)
                        {
                            if(stateControllers_[sC]->writeInTimeDir())
                            {
                                // directory: case/<timeDir>/uniform/controllers/poly/<stateControllerModel>
                                fileName sCTimePath(polyStateControllersTimePath/sCNames_[sC]);
    
                                if(!isDir(sCTimePath))
                                {
                                    mkDir(sCTimePath);
                                }
    
                                //- creating directory for different zones but of the same model
                                const word& regionName = stateControllers_[sC]->regionName();
    
                                // directory: case/<timeDir>/uniform/controllers/poly/<stateControllerModel>/<cellZoneName>
                                fileName zoneTimePath(sCTimePath/regionName);
    
                                if (!isDir(zoneTimePath))
                                {
                                    mkDir(zoneTimePath);
                                }
    
                                timePathNames[sC] = zoneTimePath;
                            }
                        }
                    }
                }
            }

            // -- write out data (do not comment this out)
            forAll(stateControllers_, sC)
            {
                stateControllers_[sC]->output(sCFixedPathNames_[sC], timePathNames[sC]);
            }
        }

        {
            List<fileName> timePathNames(fCFixedPathNames_.size());
    
            if(nFluxControllers_ > 0)
            {
                if(Pstream::master())
                {
    
                    // directory: case/<timeDir>/uniform
                    fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
                
                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }
        
                    if(fluxControllers_.size() > 0)
                    {
                    // directory: case/<timeDir>/uniform/controllers
                        fileName controllersTimePath(uniformTimePath/"controllers");
    
                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }
    
                        // directory: case/<timeDir>/uniform/controllers/poly
                        fileName polyTimePath(controllersTimePath/"poly");
                    
                        if (!isDir(polyTimePath))
                        {
                            mkDir(polyTimePath);    
                        }
    
                        // directory: case/<timeDir>/uniform/fluxControllers
                        fileName polyControllersTimePath(polyTimePath/"fluxControllers");
                    
                        if (!isDir(polyControllersTimePath))
                        {
                            mkDir(polyControllersTimePath);    
                        }
    
                        forAll(fluxControllers_, fC)
                        {
                            if(stateControllers_[fC]->writeInTimeDir())
                            {
                                // directory: case/<timeDir>/uniform/controllers/poly/<fluxControllerModel>
                                fileName fCTimePath(polyControllersTimePath/fCNames_[fC]);
            
                                if(!isDir(fCTimePath))
                                {
                                    mkDir(fCTimePath);
                                }
            
                                const word& regionName = fluxControllers_[fC]->regionName();
    
                                // directory: case/<timeDir>/uniform/controllers/poly/<fluxControllerModel>  <faceZoneName>      
                                fileName zoneTimePath(fCTimePath/regionName);
                
                                if (!isDir(zoneTimePath))
                                {
                                    mkDir(zoneTimePath);    
                                }
    
                                timePathNames[fC] = zoneTimePath;
                            }
                        }
                    }
                }
            }
    
            // -- write out data (do not comment this out)
            forAll(fluxControllers_, fC)
            {
                fluxControllers_[fC]->output(fCFixedPathNames_[fC], timePathNames[fC]);
            }
        }

        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        IOdictionary polyControllersDict
        (
            IOobject
            (
                "controllersDict",
                time_.system(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        {
            stateControllersList_.clear();
        
            stateControllersList_ = polyControllersDict.lookup("polyStateControllers");

//             Info << "updating" << endl;        
            forAll(stateControllers_, sC)
            {
                const entry& polyControllersI = stateControllersList_[sC];
                const dictionary& polyControllersIDict = polyControllersI.dict();
//                 Info << 
                stateControllers_[sC]->updateProperties(polyControllersIDict);
            }
        }

        {
            fluxControllersList_.clear();
        
            fluxControllersList_ = polyControllersDict.lookup("polyFluxControllers");
        
            forAll(fluxControllers_, fC)
            {
                const entry& polyControllersI = fluxControllersList_[fC];
                const dictionary& polyControllersIDict = polyControllersI.dict();
    
                fluxControllers_[fC]->updateProperties(polyControllersIDict);
            }
        }
    }
}

const label& polyControllers::nStateControllers() const
{
    return nStateControllers_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
