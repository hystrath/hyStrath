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

#include "pdBoundaries.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Null Constructor
pdBoundaries::pdBoundaries
(
    Time& t,
    const polyMesh& mesh
)
:
    time_(t),
    pdBoundariesDict_
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
{
}


//- Constructor for pdFoam
pdBoundaries::pdBoundaries
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud
)
:
    time_(t),
    pdBoundariesDict_
    (
        IOobject
        (
            "boundariesDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(pdBoundariesDict_.lookup("pdPatchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(pdBoundariesDict_.lookup("pdCyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(pdBoundariesDict_.lookup("pdGeneralBoundaries")),
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

            patchBoundaryModels_[p] = autoPtr<pdPatchBoundary>
            (
                pdPatchBoundary::New(t, mesh, cloud, boundaryIDict)
            );

            patchBoundaryNames_[p] = patchBoundaryModels_[p]->type();
            patchBoundaryIds_[p] = p;
            nPatchBoundaryModels_++;
        }
    }

    checkPatchBoundaryModels(mesh);

    //- cyclic boundaries

    if(cyclicBoundaryModels_.size() > 0 )
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();

            cyclicBoundaryModels_[c] = autoPtr<pdCyclicBoundary>
            (
                pdCyclicBoundary::New(t, mesh, cloud, boundaryIDict)
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

            generalBoundaryModels_[g] = autoPtr<pdGeneralBoundary>
            (
                pdGeneralBoundary::New(t, mesh, cloud, boundaryIDict)
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

        // directory: case/boundaries/pd
        fileName pdBoundariesPath(boundariesPath/"pd");

        if( !isDir(pdBoundariesPath) )
        {
            mkDir(pdBoundariesPath);
        }

        // directory: case/boundaries/pd/patchBoundaryModels
        fileName patchBoundaryModelsPath(pdBoundariesPath/"patchBoundaryModels");

        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);
        }

        forAll(patchBoundaryModels_, p)
        {
            if(patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/pd/patchBoundaryModels/<patchBoundaryModel>
                fileName patchBoundaryModelPath(patchBoundaryModelsPath/patchBoundaryNames_[p]);

                if (!isDir(patchBoundaryModelPath))
                {
                    mkDir(patchBoundaryModelPath);
                }

                const word& patchName = patchBoundaryModels_[p]->patchName();

                // directory: case/controllers/pd/patchBoundaryModels/<patchBoundaryModel>/<patchName>
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

        // directory: case/boundaries/pd
        fileName pdBoundariesPath(boundariesPath/"pd");

        if( !isDir(pdBoundariesPath) )
        {
            mkDir(pdBoundariesPath);
        }

        // directory: case/boundaries/pd/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath(pdBoundariesPath/"cyclicBoundaryModels");

        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if(cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/pd/cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName cyclicBoundaryModelPath(cyclicBoundaryModelsPath/cyclicBoundaryNames_[c]);

                if (!isDir(cyclicBoundaryModelPath))
                {
                    mkDir(cyclicBoundaryModelPath);
                }

                const word& patchName = cyclicBoundaryModels_[c]->patchName();

                // directory: case/controllers/pd/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>
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

        // directory: case/boundaries/pd
        fileName pdBoundariesPath(boundariesPath/"pd");

        if( !isDir(pdBoundariesPath) )
        {
            mkDir(pdBoundariesPath);
        }

        // directory: case/boundaries/pd/cyclicBoundaryModels
        fileName generalBoundaryModelsPath(pdBoundariesPath/"generalBoundaryModels");

        if (!isDir(generalBoundaryModelsPath))
        {
            mkDir(generalBoundaryModelsPath);
        }

        forAll(generalBoundaryModels_, g)
        {
            if(generalBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/pd/generalBoundaryModels/<generalBoundaryModel>
                fileName generalBoundaryModelPath(generalBoundaryModelsPath/generalBoundaryNames_[g]);

                if (!isDir(generalBoundaryModelPath))
                {
                    mkDir(generalBoundaryModelPath);
                }

                const word& patchName = generalBoundaryModels_[g]->patchName();

                // directory: case/controllers/pd/generalBoundaryModels/<generalBoundaryModel>/<patchName>
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

pdBoundaries::~pdBoundaries()
{
    //Info << "pdBoundaries Destructor" << endl;
}


void pdBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
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
        FatalErrorIn("pdBoundaries::checkBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of cyclic boundary models = "  << nCyclicBoundaryModels_
            << " chosen in the boundaryiesDict are inconsistent."
            << abort(FatalError);
    }
}


void pdBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
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
                FatalErrorIn("pdBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
                    << nl
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. No of models chosen for this patch are: "
                    << nPatches  << ", in "
                    << mesh.time().system()/"pdBoundariesDict"
                    << abort(FatalError);
            }
        }
    }

//     Pout << "patchToModelId_: " << patchToModelId_ << endl;

    if(nPolyPatches != nPatchBoundaryModels_)
    {
        FatalErrorIn("pdBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)")
            << nl
            << " Number of poly-patches = "  << nPolyPatches
            << " in blockMeshDict, are not equal to the number of patch models = "
            << nPatchBoundaryModels_  << ", defined in "
            << mesh.time().system()/"pdBoundariesDict"
            << abort(FatalError);
    }
}


void pdBoundaries::updateTimeInfo()
{
    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->updateTime();
    }
}

void pdBoundaries::setInitialConfig()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->setBoundaryFields();
        patchBoundaryModels_[p]->initialConfiguration();
    }

    forAll(cyclicBoundaryModels_, c)
    {
        cyclicBoundaryModels_[c]->initialConfiguration();
    }

    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->initialConfiguration();
    }
}

void pdBoundaries::calculateProps()
{
    forAll(patchBoundaryModels_, p)
    {
        patchBoundaryModels_[p]->calculateProperties();
    }

    forAll(cyclicBoundaryModels_, c)
    {
        cyclicBoundaryModels_[c]->calculateProperties();
    }

    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->calculateProperties();
    }
}

// impose model after calculation of forces
void pdBoundaries::controlBeforeMove()
{
    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->controlParcelsBeforeMove();
    }
}


void pdBoundaries::controlBeforeCollisions()
{
    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->controlParcelsBeforeCollisions();
    }
}

void pdBoundaries::controlAfterCollisions()
{
    forAll(generalBoundaryModels_, g)
    {
        generalBoundaryModels_[g]->controlParcelsAfterCollisions();
    }
}


//- output
void pdBoundaries::outputResults()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        {
            List<fileName> timePathNames(pBFixedPathNames_.size());

            if(nPatchBoundaryModels_ > 0)
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

                // directory: case/<timeDir>/uniform/boundaries/pd
                fileName pdTimePath(boundariesTimePath/"pd");

                if (!isDir(pdTimePath))
                {
                    mkDir(pdTimePath);
                }

                // directory: case/<timeDir>/uniform/boundaries/pd/patchBoundaryModels
                fileName pdPatchBoundaryModelsTimePath(pdTimePath/"patchBoundaryModels");

                if (!isDir(pdPatchBoundaryModelsTimePath))
                {
                    mkDir(pdPatchBoundaryModelsTimePath);
                }

                forAll(patchBoundaryModels_, p)
                {
                    if
                    (
                        patchBoundaryModels_[p]->writeInTimeDir() ||
                        patchBoundaryModels_[p]->writeInCase()
                    )
                    {

//                         // directory: case/<timeDir>/uniform/controllers/pd/patchBoundaryModels/<patchBoundaryModel>
//                         fileName pBTimePath(pdPatchBoundaryModelsTimePath/pBFixedPathNames_[p]);
//
//                         if(!isDir(pBTimePath))
//                         {
//                             mkDir(pBTimePath);
//                         }

                        //- creating directory for different zones but of the same model
                        const word& patchName = patchBoundaryModels_[p]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/pd/patchBoundaryModels/<patchBoundaryModel>/<patchName>
                        fileName patchTimePath(pdPatchBoundaryModelsTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[p] = patchTimePath;

                        patchBoundaryModels_[p]->output(pBFixedPathNames_[p], timePathNames[p]);
// 						patchBoundaryModels_[p]->writeField();
                    }
                }
            }
        }

        {
            List<fileName> timePathNames(cMFixedPathNames_.size());

            if(nCyclicBoundaryModels_ > 0)
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

                // directory: case/<timeDir>/uniform/boundaries/pd
                fileName pdTimePath(boundariesTimePath/"pd");

                if (!isDir(pdTimePath))
                {
                    mkDir(pdTimePath);
                }

                // directory: case/<timeDir>/uniform/boundaries/pd/cyclicBoundaryModels
                fileName pdCyclicBoundaryModelsTimePath(pdTimePath/"cyclicBoundaryModels");

                if (!isDir(pdCyclicBoundaryModelsTimePath))
                {
                    mkDir(pdCyclicBoundaryModelsTimePath);
                }

                forAll(cyclicBoundaryModels_, c)
                {
                    if
                    (
                        cyclicBoundaryModels_[c]->writeInTimeDir() ||
                        cyclicBoundaryModels_[c]->writeInCase()
                    )
                    {
                        // directory: case/<timeDir>/uniform/controllers/pd/cyclicBoundaryModels/<cyclicBoundaryModel>
                        fileName cMTimePath(pdCyclicBoundaryModelsTimePath/cMFixedPathNames_[c]);

                        if(!isDir(cMTimePath))
                        {
                            mkDir(cMTimePath);
                        }

                        //- creating directory for different zones but of the same model
                        const word& patchName = cyclicBoundaryModels_[c]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/pd/cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>
                        fileName patchTimePath(cMTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[c] = patchTimePath;

                        cyclicBoundaryModels_[c]->output(cMFixedPathNames_[c], timePathNames[c]);
                    }
                }
            }
        }

        {
            List<fileName> timePathNames(gMFixedPathNames_.size());

            if(nGeneralBoundaryModels_ > 0)
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

                // directory: case/<timeDir>/uniform/boundaries/pd
                fileName pdTimePath(boundariesTimePath/"pd");

                if (!isDir(pdTimePath))
                {
                    mkDir(pdTimePath);
                }

                // directory: case/<timeDir>/uniform/boundaries/pd/generalBoundaryModels
                fileName pdGeneralBoundaryModelsTimePath(pdTimePath/"generalBoundaryModels");

                if (!isDir(pdGeneralBoundaryModelsTimePath))
                {
                    mkDir(pdGeneralBoundaryModelsTimePath);
                }

                forAll(generalBoundaryModels_, g)
                {
                    if
                    (
                        generalBoundaryModels_[g]->writeInTimeDir() ||
                        generalBoundaryModels_[g]->writeInCase()
                    )
                    {
                        // directory: case/<timeDir>/uniform/controllers/pd/generalBoundaryModels/<generalBoundaryModel>
                        fileName gMTimePath(pdGeneralBoundaryModelsTimePath/gMFixedPathNames_[g]);

                        if(!isDir(gMTimePath))
                        {
                            mkDir(gMTimePath);
                        }

                        //- creating directory for different zones but of the same model
                        const word& patchName = generalBoundaryModels_[g]->patchName();

                        // directory: case/<timeDir>/uniform/controllers/pd/generalBoundaryModels/<generalBoundaryModel>/<patchName>
                        fileName patchTimePath(gMTimePath/patchName);

                        if (!isDir(patchTimePath))
                        {
                            mkDir(patchTimePath);
                        }

                        timePathNames[g] = patchTimePath;

                        generalBoundaryModels_[g]->output(gMFixedPathNames_[g], timePathNames[g]);
                    }
                }
            }
        }

        // RE-READ DICTIONARIES FOR MODIFIED PROPERTIES (RUN-TIME SELECTION)

        {
            patchBoundaryList_.clear();

            patchBoundaryList_ = pdBoundariesDict_.lookup("pdPatchBoundaries");

            forAll(patchBoundaryModels_, p)
            {
                const entry& boundaryI = patchBoundaryList_[p];
                const dictionary& boundaryIDict = boundaryI.dict();

                patchBoundaryModels_[p]->updateProperties(boundaryIDict);
            }
        }

        {
            cyclicBoundaryList_.clear();

            cyclicBoundaryList_ = pdBoundariesDict_.lookup("pdCyclicBoundaries");

            forAll(cyclicBoundaryModels_, c)
            {
                const entry& boundaryI = cyclicBoundaryList_[c];
                const dictionary& boundaryIDict = boundaryI.dict();

                cyclicBoundaryModels_[c]->updateProperties(boundaryIDict);
            }
        }

        {
            generalBoundaryList_.clear();

            generalBoundaryList_ = pdBoundariesDict_.lookup("pdGeneralBoundaries");

            forAll(generalBoundaryModels_, g)
            {
                const entry& boundaryI = generalBoundaryList_[g];
                const dictionary& boundaryIDict = boundaryI.dict();

                generalBoundaryModels_[g]->updateProperties(boundaryIDict);
            }
        }
    }
}

const label& pdBoundaries::nPatchBoundaryModels() const
{
    return nPatchBoundaryModels_;
}

const label& pdBoundaries::nCyclicBoundaryModels() const
{
    return nCyclicBoundaryModels_;
}

const label& pdBoundaries::nGeneralBoundaryModels() const
{
    return nGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
