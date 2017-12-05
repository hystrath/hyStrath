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

#include "dsmcCyclicBoundary.H"
#include "polyBoundaryMeshEntries.H"
#include "graph.H"
#include "processorCyclicPolyPatch.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcCyclicBoundary, 0);

defineRunTimeSelectionTable(dsmcCyclicBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCyclicBoundary::dsmcCyclicBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    time_(t),
//     rndGen_(cloud_.rndGen()),
    boundaryDict_(dict.subDict("cyclicBoundaryProperties")),

    patchName_(),
    patchId_(0),
    neighbPatchName_(),
    neighbPatchId_(0),

    faces_(),
    coupledFacesA_(),
    coupledFacesB_(),
//     cells_(),
    cellsA_(),
    cellsB_(),
    densities_(),
    velocities_(),
    temperatures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{

    const word patchName = boundaryDict_.lookup("patchName");
    patchName_ = patchName;

    //- confirm that the patch exists on the mesh
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchId_ == -1)
    {
        FatalErrorIn("dsmcCyclicBoundary::dsmcCyclicBoundary()")
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    if (!isA<cyclicPolyPatch>(patch))
    {
        FatalErrorIn("dsmcCyclicBoundary::dsmcCyclicBoundary()")
            << "Patch: " << patchName_ << " is not a cyclic boundary. " 
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    const cyclicPolyPatch& cyclicPatch =
                    refCast<const cyclicPolyPatch>(patch);

    neighbPatchName_ = cyclicPatch.neighbPatchName();

    neighbPatchId_ =  cyclicPatch.neighbPatchID();

//     Info << " patch name: " << patchName_
//          << " patch id: " << patchId_
//          << " neighbour patch name: " << neighbPatchName_
//          << " neighbour patch id: " << neighbPatchId_
//         << endl;

    getCoupledFaces(mesh_);
}

void dsmcCyclicBoundary::getCoupledFaces(const polyMesh& mesh)
{

    // in parallel, openFOAM replaces cyclic boundaries with processor boundaries
    // if processor boundaries coincide with cyclic boundaries. We therefore need
    // to do the following corrections to get the correct faces on the selected
    // cyclic patch.

    if (Pstream::parRun())
    {
        DynamicList<label> coupledFacesA(0);
        DynamicList<label> coupledFacesB(0);
        DynamicList<label> cellsA(0);
        DynamicList<label> cellsB(0);


        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        for(label i = 0; i < patch.size(); i++)
        {
            label globalFaceI = patch.start() + i;
            label cellI = patch.faceCells()[i];
            coupledFacesA.append(globalFaceI);
            cellsA.append(cellI);
        }

        const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

        for(label i = 0; i < patchN.size(); i++)
        {
            label globalFaceI = patchN.start() + i;
            label cellI = patchN.faceCells()[i];
            coupledFacesB.append(globalFaceI);
            cellsB.append(cellI);
        }

        // check for process`or-cyclic poly patches

        forAll( mesh_.boundaryMesh(), patchI )
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

//             if (isA<processorPolyPatch>(patch))
//             {
                if(isA<processorCyclicPolyPatch>(patch))
                {
                    word patchName = refCast<const processorCyclicPolyPatch>
                    (
                        patch
                    ).referPatchName();
    
                    if(patchName == patchName_)
                    {
                        for(label i = 0; i < patch.size(); i++)
                        {
                            label globalFaceI = patch.start() + i;
                            label cellI = patch.faceCells()[i];
                            coupledFacesA.append(globalFaceI);
                            cellsA.append(cellI);
                        }
                    }
                    else if(patchName == neighbPatchName_)
                    {
                        for(label i = 0; i < patch.size(); i++)
                        {
                            label globalFaceI = patch.start() + i;
                            label cellI = patch.faceCells()[i];
                            coupledFacesB.append(globalFaceI);
                            cellsB.append(cellI);
                        }
                    }
                }
//             }
        }

        coupledFacesA.shrink();
        coupledFacesB.shrink();
        cellsA.shrink();
        cellsB.shrink();

        coupledFacesA_.transfer(coupledFacesA);
        coupledFacesB_.transfer(coupledFacesB);
        cellsA_.transfer(cellsA);
        cellsB_.transfer(cellsB);
    }
    else
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        coupledFacesA_.setSize(patch.size(), -1);
        cellsA_.setSize(patch.size(), -1);

        for(label i = 0; i < patch.size(); i++)
        {
            label globalFaceI = patch.start() + i;
            label cellI = patch.faceCells()[i];
            coupledFacesA_[i] = globalFaceI;
            cellsA_[i] = cellI;
        }

        const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

        coupledFacesB_.setSize(patchN.size(), -1);
        cellsB_.setSize(patchN.size(), -1);

        for(label i = 0; i < patchN.size(); i++)
        {
            label globalFaceI = patchN.start() + i;
            label cellI = patchN.faceCells()[i];
            coupledFacesB_[i] = globalFaceI;
            cellsB_[i] = cellI;
        }
    }

    label totalFaces = coupledFacesA_.size() + coupledFacesB_.size();

    faces_.setSize(totalFaces, -1);

    label index = 0;

    forAll(coupledFacesA_, fA)
    {
        faces_[index] = coupledFacesA_[fA];
        index++;
    }

    forAll(coupledFacesB_, fB)
    {
        faces_[index] = coupledFacesB_[fB];
        index++;
    }

//     Info << "coupledFacesA: "  << coupledFacesA_ << endl;
//     Info << "coupledFacesB: "  << coupledFacesB_ << endl;
//     Info << "totalFaces: "  << faces_ << endl;

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcCyclicBoundary> dsmcCyclicBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcCyclicBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting dsmcCyclicBoundaryModel "
         << dsmcCyclicBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcCyclicBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcCyclicBoundary::New(const dictionary&) : " << endl
            << "    unknown dsmcCyclicBoundary type "
            << dsmcCyclicBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcCyclicBoundary>
	(
		cstrIter()(t, mesh, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCyclicBoundary::~dsmcCyclicBoundary()
{}

void dsmcCyclicBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("cyclicBoundaryProperties");
}

void dsmcCyclicBoundary::setNewBoundaryFields()
{
const word patchName = boundaryDict_.lookup("patchName");
    patchName_ = patchName;

    //- confirm that the patch exists on the mesh
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

//     if(patchId_ == -1)
//     {
//         FatalErrorIn("dsmcCyclicBoundary::dsmcCyclicBoundary()")
//             << "Cannot find patch: " << patchName_ << nl << "in: "
//             << t.system()/"boundariesDict"
//             << exit(FatalError);
//     }

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

//     if (!isA<cyclicPolyPatch>(patch))
//     {
//         FatalErrorIn("dsmcCyclicBoundary::dsmcCyclicBoundary()")
//             << "Patch: " << patchName_ << " is not a cyclic boundary. " 
//             << nl << "in: "
//             << t.system()/"boundariesDict"
//             << exit(FatalError);
//     }

    const cyclicPolyPatch& cyclicPatch =
                    refCast<const cyclicPolyPatch>(patch);

    neighbPatchName_ = cyclicPatch.neighbPatchName();

    neighbPatchId_ =  cyclicPatch.neighbPatchID();

//     Info << " patch name: " << patchName_
//          << " patch id: " << patchId_
//          << " neighbour patch name: " << neighbPatchName_
//          << " neighbour patch id: " << neighbPatchId_
//         << endl;

    getCoupledFaces(mesh_);
}


void dsmcCyclicBoundary::writeTimeData
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


void dsmcCyclicBoundary::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if(file.good())
    {
        forAll(yData, n)
        {
            file 
                << xData[n] << "\t" 
                << yData[n].x() << "\t" << yData[n].y() 
                << "\t" << yData[n].z() 
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void patchBoundary::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

const labelList& dsmcCyclicBoundary::controlPatch() const
{
    return coupledFacesA_;
}

const labelList& dsmcCyclicBoundary::controlZone() const
{
    return cellsA_;
}

const labelList& dsmcCyclicBoundary::allFaces() const
{
    return faces_;
}

const word& dsmcCyclicBoundary::patchName() const
{
    return patchName_;
}

const label& dsmcCyclicBoundary::patchId() const
{
    return patchId_;
}

const scalarField& dsmcCyclicBoundary::densityField() const
{
    return densities_;
}

scalarField& dsmcCyclicBoundary::densityField()
{
    return densities_;
}

const vectorField& dsmcCyclicBoundary::velocityField() const
{
    return velocities_;
}
vectorField& dsmcCyclicBoundary::velocityField()
{
    return velocities_;
}

const scalarField& dsmcCyclicBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& dsmcCyclicBoundary::temperatureField()
{
    return temperatures_;
}

const tensor& dsmcCyclicBoundary::strainRate() const
{
    return strainRate_;
}

tensor& dsmcCyclicBoundary::strainRate()
{
    return strainRate_;
}


const bool& dsmcCyclicBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcCyclicBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
