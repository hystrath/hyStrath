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

#include "cyclicBoundaries.H"
#include "cyclicPolyPatch.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Null Constructor 
cyclicBoundaries::cyclicBoundaries
(
    Time& t,
    const polyMesh& mesh
)
:    
    time_(t),
    boundariesDict_
    (
        IOobject
        (
            "cyclicBoundaries",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    cyclicBoundaryList_(),
    cyclicBoundaryNames_(),
    cyclicBoundaryIds_(),
    cMFixedPathNames_(),
    cyclicBoundaryModels_(),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1)
{}


//- General Constructor
cyclicBoundaries::cyclicBoundaries
(
    Time& t,
    const polyMesh& mesh,
    const label& dummy
)
:
    time_(t),
    boundariesDict_
    (
        IOobject
        (
            "cyclicBoundaries",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),    
    cyclicBoundaryList_(boundariesDict_.lookup("cyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1)
{
    Info << nl << "Creating the cyclic boundaries. " << nl << endl;

    //- cyclic boundaries
    const polyBoundaryMesh& bM = mesh.boundaryMesh();

    label nCyclicBoundaries = 0;

    DynamicList<word> cyclicNames(0);
    
    forAll(bM, patchI)
    {
        const polyPatch& patch = bM[patchI];

        if (isA<cyclicPolyPatch>(patch))
        {
            nCyclicBoundaries++;
            cyclicNames.append(bM.names()[patchI]);
        }
    }

    if(cyclicBoundaryModels_.size() != nCyclicBoundaries)
    {
        FatalErrorIn("cyclicBoundaries::cyclicBoundaries()")
            << "Number of cyclic boundaries defined in cyclicBoundaries dictionary: "
            << cyclicBoundaryModels_.size()
            << ", not equal to the number of cyclic boundaries in blockMesh: "
            << nCyclicBoundaries
            << nl
            << abort(FatalError);
    }

    if(cyclicBoundaryModels_.size() > 0 )
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();
    
            cyclicBoundaryModels_[c] = autoPtr<cyclicBoundary>
            (
                cyclicBoundary::New(t, mesh, boundaryIDict)
            );
    
            cyclicBoundaryNames_[c] = cyclicBoundaryModels_[c]->type();
            cyclicBoundaryIds_[c] = c;
        }
    }
    
    forAll(cyclicBoundaryModels_, i)
    {
        
        forAll(cyclicBoundaryModels_, j)
        {
            if(i != j)
            {
                bool neighbour =
                    cyclicBoundaryModels_[i]->isPatchNeighbour
                    (
                        cyclicBoundaryModels_[j]->patchName()
                    );
                                
                if(neighbour)
                {
                    cyclicBoundaryModels_[i]->setCoupledPatchInfo
                    (
                        cyclicBoundaryModels_[j]->centroid(),
                        cyclicBoundaryModels_[j]->normal()
                    );
                }
            }
        }
    }
}

cyclicBoundaries::~cyclicBoundaries()
{}




// const label& cyclicBoundaries::nCyclicBoundaryModels() const
// {
//     return nCyclicBoundaryModels_;
// }



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
