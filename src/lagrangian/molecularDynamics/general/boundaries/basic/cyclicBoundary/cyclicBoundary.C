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

#include "cyclicBoundary.H"
#include "polyBoundaryMeshEntries.H"
#include "cyclicPolyPatch.H"
// #include "graph.H"
#include "boundedBox.H" 

namespace Foam
{

    // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicBoundary, 0);

defineRunTimeSelectionTable(cyclicBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cyclicBoundary::cyclicBoundary
(
    Time& t,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    boundaryDict_(dict.subDict("cyclicBoundaryProperties")),
    patchName_(),
    patchId_(0),
//     separation_(0.0),
    separationVector_(vector::zero),
    
    theta_(constant::mathematical::pi),
    rotate_(false)
//     rotationAxis_(vector::zero),
//     rotationPt_(vector::zero),
//     RAB_(tensor::zero),
//     RBA_(tensor::zero),
{
    const word patchName = boundaryDict_.lookup("patchName");
    patchName_ = patchName;

    Info<< nl << "constructing cyclic boundary: "
        << patchName_ << nl << endl;

    //- confirm that the patch exists on the mesh
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchId_ == -1)
    {
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Cannot find patch on mesh: " << patchName_ 
            << exit(FatalError);
    }
    
    /*
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];


    if (!isA<cyclicPolyPatch>(patch))
    {
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Patch: " << patchName_ << " is not a cyclic boundary. " 
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }
 
    if (isA<coupledPolyPatch>(patch))
    {
        const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(patch);
        
        if (cpp.separated())
        {
            separationVector_ = cpp.separation()[0];
        }
    }*/

    
    boundaryPoints_ = List<vector>(boundaryDict_.lookup("boundaryPoints"));    

    if(boundaryPoints_.size() < 4)
    {
        FatalErrorIn ("cyclicBoundary.C")
            << nl << "boundaryPointsA has to be greater than 3 points. " << nl
            << abort(FatalError);
    }

    centroid_ = vector::zero;
    
    forAll(boundaryPoints_, p)
    {
        centroid_ += boundaryPoints_[p];
    }

    centroid_ /= boundaryPoints_.size();
    
    nF_ = boundaryDict_.lookup("normal");
    
    nFaces_ = (readLabel(boundaryDict_.lookup("nFaces")));
    
    
    //- Neighbour patch
    
    const word neighbPatchName = boundaryDict_.lookup("neighbourPatchName");
    
    patchNameN_ = neighbPatchName;
    
    patchNId_ = mesh_.boundaryMesh().findPatchID(patchNameN_);    

    if(patchNId_ == -1)
    {
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Cannot find neighbouring patch on mesh: " << patchNameN_ 
            << exit(FatalError);
    }
    
    Info << "neighbour patch name = " << patchNameN_ << endl;
    Info << "neighbour patch id = " << patchNId_ << endl;
    
//     const polyPatch& patchN = mesh_.boundaryMesh()[patchNId_];


}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<cyclicBoundary> cyclicBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    word cyclicBoundaryName
    (
        dict.lookup("cyclicBoundaryModel")
    );

    Info<< "Selecting cyclicBoundaryModel "
         << cyclicBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(cyclicBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "cyclicBoundary::New(const dictionary&) : " << endl
            << "    unknown cyclicBoundaryModel type "
            << cyclicBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<cyclicBoundary>
	(
		cstrIter()(t, mesh, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicBoundary::~cyclicBoundary()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cyclicBoundary::setCoupledPatchInfo
(
    const vector& centroid,
    const vector& nF
)
{
    separationVector_ = centroid - centroid_;
    
//     Info << "separationVector CPP = " << separationVector_
//          << ", new separationVector = " << separationVector 
//          << endl;
         
    theta_ = acos(nF_ & nF); 
//     Info << "theta_ " << theta_ << endl;
    if
    (
           (theta_ < (constant::mathematical::pi - SMALL))
        || (theta_ > (constant::mathematical::pi + SMALL))
    )
    {
        rotate_ = true;
        
        FatalErrorIn("cyclicBoundary::cyclicBoundary()")
            << "Patch: " << patchName_ << "is a rotational cyclic. " 
            << " Angle = " << theta_
            << "Rotational cyclics have not been implemented" 
            << exit(FatalError);
    }
}






const word& cyclicBoundary::patchName() const
{
    return patchName_;
}

const vector& cyclicBoundary::separationVector() const
{
    return separationVector_;
}

const label& cyclicBoundary::patchId() const
{
    return patchId_;
}

const List<vector>& cyclicBoundary::boundaryPoints() const
{
    return boundaryPoints_;
}

const vector& cyclicBoundary::normal() const
{
    return nF_;
}
const vector& cyclicBoundary::centroid() const
{
    return centroid_;
}

const scalar& cyclicBoundary::theta() const
{
    return theta_;
}
const bool& cyclicBoundary::rotate() const
{
    return rotate_;
}

bool cyclicBoundary::isPatchNeighbour(const word& patchNameNeighbour)
{
    if(patchNameNeighbour == patchNameN_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

const word& cyclicBoundary::patchNameN() const
{
    return patchNameN_;
}

const label& cyclicBoundary::patchIdN() const 
{
    return patchNId_;
}

/*
const vector& cyclicBoundary::normalA() const
{
    return nA_;
}
const vector& cyclicBoundary::normalB() const
{
    return nB_;
}

const label& cyclicBoundary::nFaces() const
{
    return nFaces_;
}

const labelList& cyclicBoundary::controlPatch() const
{
    return faces_;
}

const labelList& cyclicBoundary::controlPatchA() const
{
    return coupledFacesA_;
}

const labelList& cyclicBoundary::controlPatchB() const
{
    return coupledFacesB_;
}

const labelList& cyclicBoundary::controlZone() const
{
    return cells_;
}

const labelList& cyclicBoundary::controlZoneA() const
{
    return cellsA_;
}

const labelList& cyclicBoundary::controlZoneB() const
{
    return cellsB_;
}*/


/*


const vector& cyclicBoundary::cyclicTranslationVector() const
{
    return cyclicTranslationVector_;
}

const vector& cyclicBoundary::rotationAxis() const
{
    return rotationAxis_;
}

const vector& cyclicBoundary::rotationAxisStartPt() const
{
    return rotationStartPoint_;
}

const vector& cyclicBoundary::rotationPt() const
{
    return rotationPt_;
}
const tensor& cyclicBoundary::RAB() const
{
    return RAB_;
}
const tensor& cyclicBoundary::RBA() const
{
    return RBA_;
}*/



} // End namespace Foam

// ************************************************************************* //
