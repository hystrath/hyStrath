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

#include "polyGeneralBoundary.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyGeneralBoundary, 0);

defineRunTimeSelectionTable(polyGeneralBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyGeneralBoundary::polyGeneralBoundary
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    molCloud_(molCloud),
    boundaryDict_(dict.subDict("generalBoundaryProperties")),
	time_(t),
    patchName_(boundaryDict_.lookup("patchName")),
    patchId_(0),
    faces_(),
    nFaces_(0),
    cells_(),
    densities_(),
    velocities_(),
    temperatures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    //- confirm that the patch exists on the mesh

    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchId_ == -1)
    {
        FatalErrorIn("polyPatchBoundary::polyPatchBoundary()")
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

//     Pout << "patch name: " << patchName_ << ", patch size: " << patch.size() << endl;

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
    }

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyGeneralBoundary> polyGeneralBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyGeneralBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting polyGeneralBoundaryModel "
         << polyGeneralBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyGeneralBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyGeneralBoundary::New(const dictionary&) : " << endl
            << "    unknown polyGeneralBoundary type "
            << polyGeneralBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyGeneralBoundary>
	(
		cstrIter()(t, mesh, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyGeneralBoundary::~polyGeneralBoundary()
{}


void polyGeneralBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("generalBoundaryProperties");
}

const labelList& polyGeneralBoundary::controlPatch() const
{
    return faces_;
}

const labelList& polyGeneralBoundary::controlZone() const
{
    return cells_;
}

const word& polyGeneralBoundary::patchName() const
{
    return patchName_;
}

const label& polyGeneralBoundary::patchId() const
{
    return patchId_;
}


const scalarField& polyGeneralBoundary::densityField() const
{
    return densities_;
}

scalarField& polyGeneralBoundary::densityField()
{
    return densities_;
}

const vectorField& polyGeneralBoundary::velocityField() const
{
    return velocities_;
}
vectorField& polyGeneralBoundary::velocityField()
{
    return velocities_;
}

const scalarField& polyGeneralBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& polyGeneralBoundary::temperatureField()
{
    return temperatures_;
}

const bool& polyGeneralBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyGeneralBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
