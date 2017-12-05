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

#include "polyCyclicBoundary.H"
#include "polyBoundaryMeshEntries.H"
#include "graph.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyCyclicBoundary, 0);

defineRunTimeSelectionTable(polyCyclicBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyCyclicBoundary::polyCyclicBoundary
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    molCloud_(molCloud),
    time_(t),
    patchName_(),
    patchId_(0),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const word patchName = dict.lookup("patchName");
    patchName_ = patchName;

    Info << nl << "constructing cyclic boundary: " << patchName << nl << endl;

    //- confirm that the patch exists on the mesh
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchId_ == -1)
    {
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }


    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    if (!isA<cyclicPolyPatch>(patch))
    {
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Patch: " << patchName_ << " is not a cyclic boundary. " 
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyCyclicBoundary> polyCyclicBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyCyclicBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting polyCyclicBoundaryModel "
         << polyCyclicBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyCyclicBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyCyclicBoundary::New(const dictionary&) : " << endl
            << "    unknown polyCyclicBoundary type "
            << polyCyclicBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyCyclicBoundary>
	(
		cstrIter()(t, mesh, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyCyclicBoundary::~polyCyclicBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyCyclicBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{}

const word& polyCyclicBoundary::patchName() const
{
    return patchName_;
}

const label& polyCyclicBoundary::patchId() const
{
    return patchId_;
}

const bool& polyCyclicBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyCyclicBoundary::writeInCase() const
{
    return writeInCase_;
}

} // End namespace Foam

// ************************************************************************* //
