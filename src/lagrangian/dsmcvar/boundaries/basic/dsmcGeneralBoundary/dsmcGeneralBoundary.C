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

#include "dsmcGeneralBoundary.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcGeneralBoundary, 0);

defineRunTimeSelectionTable(dsmcGeneralBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcGeneralBoundary::dsmcGeneralBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    boundaryDict_(dict.subDict("generalBoundaryProperties")),
//     timeDict_(boundaryDict_.subDict("timeProperties")),
//     time_(t, timeDict_),
    patchName_(boundaryDict_.lookup("patchName")),
    patchId_(0),
    faces_(),
    nFaces_(0),
    patchSurfaceArea_(0.0),
    cells_(),
    density_(0.0),
    velocity_(vector::zero),
    temperature_(0.0),
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
        FatalErrorIn("atomisticPatchBoundary::atomisticPatchBoundary()")
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
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcGeneralBoundary> dsmcGeneralBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcGeneralBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting dsmcGeneralBoundaryModel "
         << dsmcGeneralBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcGeneralBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcGeneralBoundary::New(const dictionary&) : " << endl
            << "    unknown dsmcGeneralBoundary type "
            << dsmcGeneralBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<dsmcGeneralBoundary>
    (
        cstrIter()(t, mesh, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcGeneralBoundary::~dsmcGeneralBoundary()
{}

void dsmcGeneralBoundary::updateTime()
{
//     time_++;
//     timeMeas_++;
// 
//     const scalar& t = time_.time().timeOutputValue();
// 
//     if((t - initialTime_) < timePeriod_)
//     {
//         time_.controlTimeInterval().endTime() = false;
// //         time_.calcPropTimeInterval().endTime() = false;
//     }
}

void dsmcGeneralBoundary::updateBoundaryProperties(const dictionary&)
{}

void dsmcGeneralBoundary::setNewBoundaryFields()
{
    //- confirm that the patch exists on the mesh

    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique
    
    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}

void dsmcGeneralBoundary::writeTimeData
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


const labelList& dsmcGeneralBoundary::controlPatch() const
{
    return faces_;
}

const labelList& dsmcGeneralBoundary::controlZone() const
{
    return cells_;
}

const word& dsmcGeneralBoundary::patchName() const
{
    return patchName_;
}

const label& dsmcGeneralBoundary::patchId() const
{
    return patchId_;
}



const scalar& dsmcGeneralBoundary::density() const
{
    return density_;
}

scalar& dsmcGeneralBoundary::density()
{
    return density_;
}

const vector& dsmcGeneralBoundary::velocity() const
{
    return velocity_;
}
vector& dsmcGeneralBoundary::velocity()
{
    return velocity_;
}

const scalar& dsmcGeneralBoundary::temperature() const
{
    return temperature_;
}

scalar& dsmcGeneralBoundary::temperature()
{
    return temperature_;
}


const scalarField& dsmcGeneralBoundary::densityField() const
{
    return densities_;
}

scalarField& dsmcGeneralBoundary::densityField()
{
    return densities_;
}

const vectorField& dsmcGeneralBoundary::velocityField() const
{
    return velocities_;
}
vectorField& dsmcGeneralBoundary::velocityField()
{
    return velocities_;
}

const scalarField& dsmcGeneralBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& dsmcGeneralBoundary::temperatureField()
{
    return temperatures_;
}

const bool& dsmcGeneralBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcGeneralBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
