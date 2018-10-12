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

#include "pdGeneralBoundary.H"
#include "graph.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdGeneralBoundary, 0);

defineRunTimeSelectionTable(pdGeneralBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdGeneralBoundary::pdGeneralBoundary
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
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
    /**************************************************************/
    //charge_(0.0)
    /**************************************************************/
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

autoPtr<pdGeneralBoundary> pdGeneralBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
{
    word pdGeneralBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting pdGeneralBoundaryModel "
         << pdGeneralBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdGeneralBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdGeneralBoundary::New(const dictionary&) : " << endl
            << "    unknown pdGeneralBoundary type "
            << pdGeneralBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc() << abort(FatalError);
    }

    return autoPtr<pdGeneralBoundary>
    (
        cstrIter()(t, mesh, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdGeneralBoundary::~pdGeneralBoundary()
{}

void pdGeneralBoundary::updateTime()
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

void pdGeneralBoundary::updateBoundaryProperties(const dictionary&)
{}


void pdGeneralBoundary::writeTimeData
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


const labelList& pdGeneralBoundary::controlPatch() const
{
    return faces_;
}

const labelList& pdGeneralBoundary::controlZone() const
{
    return cells_;
}

const word& pdGeneralBoundary::patchName() const
{
    return patchName_;
}

const label& pdGeneralBoundary::patchId() const
{
    return patchId_;
}



const scalar& pdGeneralBoundary::density() const
{
    return density_;
}

scalar& pdGeneralBoundary::density()
{
    return density_;
}

const vector& pdGeneralBoundary::velocity() const
{
    return velocity_;
}
vector& pdGeneralBoundary::velocity()
{
    return velocity_;
}

const scalar& pdGeneralBoundary::temperature() const
{
    return temperature_;
}

scalar& pdGeneralBoundary::temperature()
{
    return temperature_;
}

/**************************************************************/
/*
const scalar& pdGeneralBoundary::charge() const
{
    return charge_;
}

scalar& pdGeneralBoundary::charge()
{
    return charge_;
}
*/
/**************************************************************/

const scalarField& pdGeneralBoundary::densityField() const
{
    return densities_;
}

scalarField& pdGeneralBoundary::densityField()
{
    return densities_;
}

const vectorField& pdGeneralBoundary::velocityField() const
{
    return velocities_;
}
vectorField& pdGeneralBoundary::velocityField()
{
    return velocities_;
}

const scalarField& pdGeneralBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& pdGeneralBoundary::temperatureField()
{
    return temperatures_;
}

/**************************************************************/
/*
const scalarField& pdGeneralBoundary::temperatureField() const
{
    return charge_;
}

scalarField& pdGeneralBoundary::temperatureField()
{
    return charge_;
}
*/
/**************************************************************/

const bool& pdGeneralBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& pdGeneralBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
