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

#include "dsmcControllers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//- Null Constructor 
dsmcControllers::dsmcControllers
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    dsmcControllersDict_
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
    fluxControllers_()
{}

//- Constructor for mdFOAM
dsmcControllers::dsmcControllers
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    dsmcControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    nFluxControllers_(0),

	stateControllersList_(dsmcControllersDict_.lookup("dsmcStateControllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
	stateControllers_(stateControllersList_.size()),

    fluxControllersList_(dsmcControllersDict_.lookup("dsmcFluxControllers")),
    fCNames_(fluxControllersList_.size()),
    fCIds_(fluxControllersList_.size()),
    fCFixedPathNames_(fluxControllersList_.size()),
	fluxControllers_(fluxControllersList_.size())
{

    Info << nl << "Creating dsmcControllers" << nl << endl;

    //- state dsmcControllers

    if(stateControllers_.size() > 0 )
    {
        forAll(stateControllers_, sC)
        {
            const entry& dsmcControllersI = stateControllersList_[sC];
            const dictionary& dsmcControllersIDict = dsmcControllersI.dict();
    
            stateControllers_[sC] = autoPtr<dsmcStateController>
            (
                dsmcStateController::New(time_, cloud, dsmcControllersIDict)
            );
    
            sCNames_[sC] = stateControllers_[sC]->type();
            sCIds_[sC] = sC;
    
            nStateControllers_++;
        }
    }

    //- flux dsmcControllers

    if(fluxControllers_.size() > 0 )
    {
        forAll(fluxControllers_, fC)
        {
            const entry& dsmcControllersI = fluxControllersList_[fC];
    
            const dictionary& dsmcControllersIDict = dsmcControllersI.dict();
    
            fluxControllers_[fC] = autoPtr<dsmcFluxController>
            (
                dsmcFluxController::New(time_, cloud, dsmcControllersIDict)
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

        // directory: case/controllers/dsmc
        fileName dsmcControllersPath(controllersPath/"dsmc");

        if( !isDir(dsmcControllersPath) )
        {
            mkDir(dsmcControllersPath);
        }

        // directory: case/controllers/dsmc/stateControllers
        fileName stateControllersPath(dsmcControllersPath/"stateControllers");
    
        if (!isDir(stateControllersPath))
        {
            mkDir(stateControllersPath);    
        }

        forAll(stateControllers_, sC)
        {
            if(stateControllers_[sC]->writeInCase())
            {
                // directory: case/controllers/dsmc/stateControllers/<stateControllerModel>
                fileName stateControllerPath(stateControllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);    
                }
    
                const word& regionName = stateControllers_[sC]->regionName();

                // directory: case/controllers/dsmc/stateControllers/<stateControllerModel>/<cellZoneName>    
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

        // directory: case/controllers/dsmc
        fileName dsmcControllersPath(time_.path()/"dsmc");

        if( !isDir(dsmcControllersPath) )
        {
            mkDir(dsmcControllersPath);
        }

        // directory: case/controllers/dsmc/fluxControllers
        fileName fluxControllersPath(dsmcControllersPath/"fluxControllers");
    
        if (!isDir(fluxControllersPath))
        {
            mkDir(fluxControllersPath);    
        }

        forAll(fluxControllers_, fC)
        {
            if(fluxControllers_[fC]->writeInCase())
            {
                // directory: case/controllers/dsmc/fluxControllers/<fluxControllerModel>
                fileName fluxControllerPath(fluxControllersPath/fCNames_[fC]);
    
                if (!isDir(fluxControllerPath))
                {
                    mkDir(fluxControllerPath);    
                }

                const word& regionName = fluxControllers_[fC]->regionName();
    
                // directory: case/controllers/dsmc/fluxControllers/<fluxControllerModel>/<faceZoneName>
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

dsmcControllers::~dsmcControllers()
{}

//- initial configuration
//- call this function after the dsmcMoleculeCloud is completely initialised
void dsmcControllers::initialConfig()
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

        //- different control stages 
void dsmcControllers::controlBeforeMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsBeforeMove();
    }
}


void dsmcControllers::controlBeforeCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsBeforeCollisions();
    }
}

void dsmcControllers::controlAfterCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsAfterCollisions();
    }
}



//- calculate properties -- call this at the end of the MD time-step.
void dsmcControllers::calculateProps()
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

//- this function is to be called at the beginning of the MD time-step. 
//  since we have placed a non-referenced time-data class in the state-controller class.
void dsmcControllers::updateTimeInfo()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->updateTime();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->updateTime();
    }
}

//- output -- call this function at the end of the MD time-step
void dsmcControllers::outputResults() 
{
    const Time& runTime = time_;

    if(runTime.outputTime())
    {
        // -- creating a set of directories in the current time directory
        {
            List<fileName> timePathNames(sCFixedPathNames_.size()); 
    
            if(nStateControllers_ > 0)
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

                    // directory: case/<timeDir>/uniform/controllers/dsmc
                    fileName dsmcTimePath(controllersTimePath/"dsmc");
                
                    if (!isDir(dsmcTimePath))
                    {
                        mkDir(dsmcTimePath);    
                    }

                    // directory: case/<timeDir>/uniform/controllers/dsmc/
                    fileName dsmcStateControllersTimePath(dsmcTimePath/"stateControllers");
                
                    if (!isDir(dsmcStateControllersTimePath))
                    {
                        mkDir(dsmcStateControllersTimePath);    
                    }

                    forAll(stateControllers_, sC)
                    {
                        if
                        (
                            stateControllers_[sC]->writeInTimeDir()
                        ) 
                        {
                            // directory: case/<timeDir>/uniform/controllers/dsmc/<stateControllerModel>
                            fileName sCTimePath(dsmcStateControllersTimePath/sCNames_[sC]);

                            if(!isDir(sCTimePath))
                            {
                                mkDir(sCTimePath);
                            }

                            //- creating directory for different zones but of the same model
                            const word& regionName = stateControllers_[sC]->regionName();

                            // directory: case/<timeDir>/uniform/controllers/dsmc/<stateControllerModel>/<cellZoneName>
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

                    // directory: case/<timeDir>/uniform/controllers/dsmc
                    fileName dsmcTimePath(controllersTimePath/"dsmc");
                
                    if (!isDir(dsmcTimePath))
                    {
                        mkDir(dsmcTimePath);    
                    }

                    // directory: case/<timeDir>/uniform/fluxControllers
                    fileName dsmcControllersTimePath(dsmcTimePath/"fluxControllers");
                
                    if (!isDir(dsmcControllersTimePath))
                    {
                        mkDir(dsmcControllersTimePath);    
                    }

                    forAll(fluxControllers_, fC)
                    {
                        if
                        (
                            fluxControllers_[fC]->writeInTimeDir()
                        )
                        {
                            // directory: case/<timeDir>/uniform/controllers/dsmc/<fluxControllerModel>
                            fileName fCTimePath(dsmcControllersTimePath/fCNames_[fC]);
        
                            if(!isDir(fCTimePath))
                            {
                                mkDir(fCTimePath);
                            }
        
                            const word& regionName = fluxControllers_[fC]->regionName();

                            // directory: case/<timeDir>/uniform/controllers/dsmc/<fluxControllerModel>  <faceZoneName>      
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

            // -- write out data (do not comment this out)
            forAll(fluxControllers_, fC)
            {
                fluxControllers_[fC]->output(fCFixedPathNames_[fC], timePathNames[fC]);
            }
        }

        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        {
            stateControllersList_.clear();
        
            stateControllersList_ = dsmcControllersDict_.lookup("dsmcStateControllers");
        
            forAll(stateControllers_, sC)
            {
                const entry& dsmcControllersI = stateControllersList_[sC];
                const dictionary& dsmcControllersIDict = dsmcControllersI.dict();
    
                stateControllers_[sC]->updateProperties(dsmcControllersIDict);
            }
        }

        {
            fluxControllersList_.clear();
        
            fluxControllersList_ = dsmcControllersDict_.lookup("dsmcFluxControllers");
        
            forAll(fluxControllers_, fC)
            {
                const entry& dsmcControllersI = fluxControllersList_[fC];
                const dictionary& dsmcControllersIDict = dsmcControllersI.dict();
    
                fluxControllers_[fC]->updateProperties(dsmcControllersIDict);
            }
        }
    }
}

const label& dsmcControllers::nStateControllers() const
{
    return nStateControllers_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
