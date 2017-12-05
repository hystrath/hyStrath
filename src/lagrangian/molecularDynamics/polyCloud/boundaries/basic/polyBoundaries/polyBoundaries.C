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

#include "polyBoundaries.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Null Constructor 
polyBoundaries::polyBoundaries
(
    Time& t,
    const polyMesh& mesh
)
:    
    time_(t),
    polyBoundariesDict_
    (
        IOobject
        (
            "boundariesDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(),
    patchBoundaryNames_(),
    patchBoundaryIds_(),
    pBFixedPathNames_(),
    patchBoundaryModels_(),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(),
    cyclicBoundaryNames_(),
    cyclicBoundaryIds_(),
    cMFixedPathNames_(),
    cyclicBoundaryModels_(),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(),
    generalBoundaryNames_(),
    generalBoundaryIds_(),
    gMFixedPathNames_(),
    generalBoundaryModels_()

{}


//- Constructor for gnemdFOAM
polyBoundaries::polyBoundaries
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    time_(t),
    polyBoundariesDict_
    (
        IOobject
        (
            "boundariesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(polyBoundariesDict_.lookup("polyPatchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(polyBoundariesDict_.lookup("polyCyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(polyBoundariesDict_.lookup("polyGeneralBoundaries")),
    generalBoundaryNames_(generalBoundaryList_.size()),
    generalBoundaryIds_(generalBoundaryList_.size()),
    gMFixedPathNames_(generalBoundaryList_.size()),
    generalBoundaryModels_(generalBoundaryList_.size())
{
    Info << "Creating the boundary models: " << nl << endl;

    //- patch boundaries

    if( patchBoundaryModels_.size() > 0 )
    {
        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            patchBoundaryModels_[p] = autoPtr<polyPatchBoundary>
            (
                polyPatchBoundary::New(t, mesh, molCloud, boundaryIDict)
            );
    
            patchBoundaryNames_[p] = patchBoundaryModels_[p]->type();
            patchBoundaryIds_[p] = p;
            nPatchBoundaryModels_++;
        }
    }

    checkPatchBoundaryModels(mesh);
    
    
    //- cyclic boundaries
    
    if( cyclicBoundaryModels_.size() > 0 )
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            cyclicBoundaryModels_[c] = autoPtr<polyCyclicBoundary>
            (
                polyCyclicBoundary::New(t, mesh, molCloud, boundaryIDict)
            );
    
            cyclicBoundaryNames_[c] = cyclicBoundaryModels_[c]->type();
            cyclicBoundaryIds_[c] = c;
            nCyclicBoundaryModels_++;
        }
    }    
    
    checkCyclicBoundaryModels(mesh);

    //- general boundaries
 
    if(generalBoundaryModels_.size() > 0 )
    {
        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            generalBoundaryModels_[g] = autoPtr<polyGeneralBoundary>
            (
                polyGeneralBoundary::New(t, mesh, molCloud, boundaryIDict)
            );
    
            generalBoundaryNames_[g] = generalBoundaryModels_[g]->type();
            generalBoundaryIds_[g] = g;
            nGeneralBoundaryModels_++;
        }
    }

    //- creating directories
    if(nPatchBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/poly
        fileName polyBoundariesPath(boundariesPath/"poly");

        if( !isDir(polyBoundariesPath) )
        {
            mkDir(polyBoundariesPath);
        }

        // directory: case/boundaries/poly/patchBoundaryModels
        fileName patchBoundaryModelsPath(polyBoundariesPath/"patchBoundaryModels");
    
        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);    
        }

        forAll(patchBoundaryModels_, p)
        {
            if(patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/poly/patchBoundaryModels/<patchBoundaryModel>
                fileName patchBoundaryModelPath(patchBoundaryModelsPath/patchBoundaryNames_[p]);

                if (!isDir(patchBoundaryModelPath))
                {
                    mkDir(patchBoundaryModelPath);    
                }
    
                const word& patchName = patchBoundaryModels_[p]->patchName();

                // directory: case/controllers/poly/patchBoundaryModels/<patchBoundaryModel>/<patchName>    
                fileName patchPath(patchBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                pBFixedPathNames_[p] = patchPath;
            }
        }
    }

    //- creating directories
    if(nCyclicBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/poly
        fileName polyBoundariesPath(boundariesPath/"poly");

        if( !isDir(polyBoundariesPath) )
        {
            mkDir(polyBoundariesPath);
        }

        // directory: case/boundaries/poly/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath(polyBoundariesPath/"cyclicBoundaryModels");
    
        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);    
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if(cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/poly/cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName cyclicBoundaryModelPath(cyclicBoundaryModelsPath/cyclicBoundaryNames_[c]);

                if (!isDir(cyclicBoundaryModelPath))
                {
                    mkDir(cyclicBoundaryModelPath);    
                }
                
                const word& patchName = cyclicBoundaryModels_[c]->patchName();    

                // directory: case/controllers/poly/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>      
                fileName patchPath(cyclicBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                cMFixedPathNames_[c] = patchPath;
            }
        }
    }

    //- creating directories
    if(nGeneralBoundaryModels_ > 0) 
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if( !isDir(boundariesPath) )
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/poly
        fileName polyBoundariesPath(boundariesPath/"poly");

        if( !isDir(polyBoundariesPath) )
        {
            mkDir(polyBoundariesPath);
        }

        // directory: case/boundaries/poly/cyclicBoundaryModels
        fileName generalBoundaryModelsPath(polyBoundariesPath/"generalBoundaryModels");
    
        if (!isDir(generalBoundaryModelsPath))
        {
            mkDir(generalBoundaryModelsPath);    
        }

        forAll(generalBoundaryModels_, g)
        {
            if(generalBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/poly/generalBoundaryModels/<generalBoundaryModel>
                fileName generalBoundaryModelPath(generalBoundaryModelsPath/generalBoundaryNames_[g]);

                if (!isDir(generalBoundaryModelPath))
                {
                    mkDir(generalBoundaryModelPath);    
                }
    
                const word& patchName = generalBoundaryModels_[g]->patchName();

                // directory: case/controllers/poly/generalBoundaryModels/<generalBoundaryModel>/<patchName>      
                fileName patchPath(generalBoundaryModelPath/patchName);
   
                if (!isDir(patchPath))
                {
                    mkDir(patchPath);    
                }
    
                gMFixedPathNames_[g] = patchPath;
            }
        }
    }

}

polyBoundaries::~polyBoundaries()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
{
    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if(isA<cyclicPolyPatch>(patch))
        {
            label patchIndex = patch.index();

            forAll(cyclicBoundaryModels_, c)
            {
                const label& patchId = cyclicBoundaryModels_[c]->patchId();
 
                if(patchIndex == patchId)
                {
                    nPolyPatches++;
                    cyclicBoundaryToModelId_[patchi] = c;
                }
            }
        }
    }

    if(nPolyPatches != nCyclicBoundaryModels_)
    {
        FatalErrorIn("polyBoundaries::checkBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of cyclic boundary models = "  << nCyclicBoundaryModels_ 
            << " chosen in the boundaryiesDict are inconsistent." 
            << abort(FatalError);
    }
}


void polyBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
{
    //- check that all poly-patches defined within blockMeshDict,
    //  each have one model.

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
    
        if
        (
            isA<polyPatch>(patch) &&
            !isA<cyclicPolyPatch>(patch) &&
            !isA<processorPolyPatch>(patch) &&
            !isA<wallPolyPatch>(patch) &&
            !isA<emptyPolyPatch>(patch) &&
            !isA<symmetryPolyPatch>(patch)
        )
        {
            nPolyPatches++;

            label patchIndex = patch.index();

            label nPatches = 0;

            forAll(patchBoundaryModels_, p)
            {
                const label& patchId = patchBoundaryModels_[p]->patchId();
 
                if(patchIndex == patchId)
                {
                    nPatches++;
                    patchToModelId_[patchi] = p;
                }
            }

            if(nPatches > 1)
            {
                FatalErrorIn("polyBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
                    << nl
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. No of models chosen for this patch are: " 
                    << nPatches  << ", in " 
                    << mesh.time().system()/"polyBoundariesDict"
                    << abort(FatalError);
            }
        }
    }

//     Pout << "patchToModelId_: " << patchToModelId_ << endl;

    if(nPolyPatches != nPatchBoundaryModels_)
    {
        FatalErrorIn("polyBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of poly-patches = "  << nPolyPatches 
            << " in blockMeshDict, are not equal to the number of patch models = " 
            << nPatchBoundaryModels_  << ", defined in " 
            << mesh.time().system()/"polyBoundariesDict"
            << abort(FatalError);
    }
}




void polyBoundaries::setInitialConfig()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->initialConfiguration();
    }

//     forAll(cyclicBoundaryModels_, c)
//     {
//         cyclicBoundaryModels_[c]->initialConfiguration();
//     }

    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->initialConfiguration();
    }
}

void polyBoundaries::calculateProps()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->calculateProperties();
    }

//     forAll(cyclicBoundaryModels_, c)
//     {
//         cyclicBoundaryModels_[c]->calculateProperties();
//     }

    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->calculateProperties();
    }
}

void polyBoundaries::controlAfterMove()
{
//     forAll(cyclicBoundaryModels_, c)
//     {
//         cyclicBoundaryModels_[c]->controlAfterMove();
//     }
}

// impose model after calculation of forces
void polyBoundaries::controlAfterForces()
{
    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->controlMols();
    }
}

//- output
void polyBoundaries::outputResults()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        //- PATCH BOUNDARY MODELS
        {
            List<fileName> timePathNames(pBFixedPathNames_.size()); 
    
            if(nPatchBoundaryModels_ > 0)
            {
                if(Pstream::master())
                {
                    // directory: case/<timeDir>/uniform
                    fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
                
                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries
                    fileName boundariesTimePath(uniformTimePath/"boundaries");
    
                    if (!isDir(boundariesTimePath))
                    {
                        mkDir(boundariesTimePath);
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries/poly
                    fileName polyTimePath(boundariesTimePath/"poly");
                
                    if (!isDir(polyTimePath))
                    {
                        mkDir(polyTimePath);    
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries/poly/patchBoundaryModels
                    fileName polyPatchBoundaryModelsTimePath(polyTimePath/"patchBoundaryModels");
                
                    if (!isDir(polyPatchBoundaryModelsTimePath))
                    {
                        mkDir(polyPatchBoundaryModelsTimePath);    
                    }
    
                    forAll(patchBoundaryModels_, p)
                    {
                        if(patchBoundaryModels_[p]->writeInTimeDir())
                        {
                            // directory: case/<timeDir>/uniform/controllers/poly/patchBoundaryModels/<patchBoundaryModel>
                            fileName pBTimePath(polyPatchBoundaryModelsTimePath/pBFixedPathNames_[p]);
    
                            if(!isDir(pBTimePath))
                            {
                                mkDir(pBTimePath);
                            }
    
                            //- creating directory for different zones but of the same model
                            const word& patchName = patchBoundaryModels_[p]->patchName();
    
                            // directory: case/<timeDir>/uniform/controllers/poly/patchBoundaryModels/<patchBoundaryModel>/<patchName>
                            fileName patchTimePath(pBTimePath/patchName);
    
                            if (!isDir(patchTimePath))
                            {
                                mkDir(patchTimePath);
                            }
    
                            timePathNames[p] = patchTimePath;
                        }
                    }
                }
            }

            forAll(patchBoundaryModels_, p)
            {
                patchBoundaryModels_[p]->output(pBFixedPathNames_[p], timePathNames[p]);
            }
        }

        //- GENERAL BOUNDARY MODELS
        {
            List<fileName> timePathNames(gMFixedPathNames_.size()); 
    
            if(nGeneralBoundaryModels_ > 0)
            {
                if(Pstream::master())
                {
                    // directory: case/<timeDir>/uniform
                    fileName uniformTimePath(runTime.path()/runTime.timeName()/"uniform");
                
                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries
                    fileName boundariesTimePath(uniformTimePath/"boundaries");
    
                    if (!isDir(boundariesTimePath))
                    {
                        mkDir(boundariesTimePath);
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries/poly
                    fileName polyTimePath(boundariesTimePath/"poly");
                
                    if (!isDir(polyTimePath))
                    {
                        mkDir(polyTimePath);    
                    }
    
                    // directory: case/<timeDir>/uniform/boundaries/poly/patchBoundaryModels
                    fileName polyGeneralBoundaryModelsTimePath
                    (
                        polyTimePath/"generalBoundaryModels"
                    );
                
                    if (!isDir(polyGeneralBoundaryModelsTimePath))
                    {
                        mkDir(polyGeneralBoundaryModelsTimePath);    
                    }
    
                    forAll(generalBoundaryModels_, g)
                    {
                        if
                        (
                            generalBoundaryModels_[g]->writeInTimeDir()
                        )
                        {
                            /* directory: case/<timeDir>/uniform/controllers/poly/
                            patchBoundaryModels/<patchBoundaryModel>*/
                            fileName gBTimePath(polyGeneralBoundaryModelsTimePath/gMFixedPathNames_[g]);
    
                            if(!isDir(gBTimePath))
                            {
                                mkDir(gBTimePath);
                            }
    
                            //- creating directory for different zones but of the same model
                            const word& patchName = generalBoundaryModels_[g]->patchName();
    
                            /* directory: case/<timeDir>/uniform/controllers/poly/
                                patchBoundaryModels/<patchBoundaryModel>/<patchName> */
                            fileName patchTimePath(gBTimePath/patchName);
    
                            if (!isDir(patchTimePath))
                            {
                                mkDir(patchTimePath);
                            }
    
                            timePathNames[g] = patchTimePath;
                        }
                    }
                }
            }

            // -- write out data (do not comment this out)
            forAll(generalBoundaryModels_, g)
            {
                generalBoundaryModels_[g]->output(gMFixedPathNames_[g], timePathNames[g]);
            }
        }


        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        {
            patchBoundaryList_.clear();
        
            patchBoundaryList_ = polyBoundariesDict_.lookup("polyPatchBoundaries");
        
            forAll(patchBoundaryModels_, p)
            {
                const entry& boundaryI = patchBoundaryList_[p];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                patchBoundaryModels_[p]->updateProperties(boundaryIDict);
            }
        }
/*
        {
            cyclicBoundaryList_.clear();
        
            cyclicBoundaryList_ = polyBoundariesDict_.lookup("polyCyclicBoundaries");
        
            forAll(cyclicBoundaryModels_, c)
            {
                const entry& boundaryI = cyclicBoundaryList_[c];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                cyclicBoundaryModels_[c]->updateProperties(boundaryIDict);
            }
        }*/

        {
            generalBoundaryList_.clear();
        
            generalBoundaryList_ = polyBoundariesDict_.lookup("polyGeneralBoundaries");
        
            forAll(generalBoundaryModels_, g)
            {
                const entry& boundaryI = generalBoundaryList_[g];
                const dictionary& boundaryIDict = boundaryI.dict();
    
                generalBoundaryModels_[g]->updateProperties(boundaryIDict);
            }
        }
    }
}



const label& polyBoundaries::nPatchBoundaryModels() const
{
    return nPatchBoundaryModels_;
}

// const label& polyBoundaries::nCyclicBoundaryModels() const
// {
//     return nCyclicBoundaryModels_;
// }

const label& polyBoundaries::nGeneralBoundaryModels() const
{
    return nGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
