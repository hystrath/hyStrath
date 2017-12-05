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

#include "polyPatchBoundary.H"
#include "IFstream.H"
#include "graph.H"
#include "polyMoleculeCloud.H"
namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyPatchBoundary, 0);

defineRunTimeSelectionTable(polyPatchBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPatchBoundary::polyPatchBoundary
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    t_(t),
    molCloud_(molCloud),
    boundaryDict_(dict.subDict("patchBoundaryProperties")),
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

//     Pout << "faces: " << faces_ << endl;

    //- set the targeted fields
//     densities_.setSize(cells_.size(), scalar(0.0));
//     velocities_.setSize(cells_.size(), vector::zero);
//     temperatures_.setSize(cells_.size(), scalar(0.0));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<polyPatchBoundary> polyPatchBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
{
    word polyPatchBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting polyPatchBoundaryModel "
         << polyPatchBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(polyPatchBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "polyPatchBoundary::New(const dictionary&) : " << endl
            << "    unknown polyPatchBoundary type "
            << polyPatchBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<polyPatchBoundary>
	(
		cstrIter()(t, mesh, molCloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPatchBoundary::~polyPatchBoundary()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- needs to be called after the referred interaction list is built in polyMoleculeCloud
void polyPatchBoundary::setBoundaryCells()
{}

void polyPatchBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("patchBoundaryProperties");
}

const labelList& polyPatchBoundary::controlPatch() const
{
    return faces_;
}

const labelList& polyPatchBoundary::controlZone() const
{
    return cells_;
}

const word& polyPatchBoundary::patchName() const
{
    return patchName_;
}

const label& polyPatchBoundary::patchId() const
{
    return patchId_;
}


const scalarField& polyPatchBoundary::densityField() const
{
    return densities_;
}

scalarField& polyPatchBoundary::densityField()
{
    return densities_;
}

const vectorField& polyPatchBoundary::velocityField() const
{
    return velocities_;
}
vectorField& polyPatchBoundary::velocityField()
{
    return velocities_;
}

const scalarField& polyPatchBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& polyPatchBoundary::temperatureField()
{
    return temperatures_;
}

const bool& polyPatchBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& polyPatchBoundary::writeInCase() const
{
    return writeInCase_;
}


scalar polyPatchBoundary::avReqDensity()
{
    scalar totalDensity = 0.0;

    forAll(densities_, c)
    {
        totalDensity += densities_[c];
    }

    if(cells_.size() > 0)
    {
        totalDensity /= scalar(cells_.size());
    }

    return totalDensity;
}

vector polyPatchBoundary::avReqVelocity()
{
    vector totalVel = vector::zero;

    forAll(velocities_, c)
    {
        totalVel += velocities_[c];
    }

    if(cells_.size() > 0)
    {
        totalVel /= scalar(cells_.size());
    }

    return totalVel;
}

scalar polyPatchBoundary::avReqTemperature()
{
    scalar totalTemp = 0.0;

    forAll(densities_, c)
    {
        totalTemp += temperatures_[c];
    }

    if(cells_.size() > 0)
    {
        totalTemp /= scalar(cells_.size());
    }

    return totalTemp;
}



} // End namespace Foam

// ************************************************************************* //
