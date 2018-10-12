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

#include "pdPatchBoundary.H"
#include "IFstream.H"
#include "graph.H"
#include "pdCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pdPatchBoundary, 0);

defineRunTimeSelectionTable(pdPatchBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdPatchBoundary::pdPatchBoundary
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
//     rndGen_(clock::getTime()),
    boundaryDict_(dict.subDict("patchBoundaryProperties")),
    patchName_(boundaryDict_.lookup("patchName")),
    patchId_(0),
    faces_(),
    nFaces_(0),
    cells_(),
    patchSurfaceArea_(0.0),
    density_(0.0),
    velocity_(vector::zero),
    temperature_(0.0),
    densities_(),
    velocities_(),
    temperatures_(),
    writeInTimeDir_(true),
    writeInCase_(true),
    measurePropertiesAtWall_(false),
    preIE_(0.0),
    preIMom_(vector::zero)
{
    //- confirm that the patch exists on the mesh

    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if(patchId_ == -1)
    {
        FatalErrorIn("pdPatchBoundary::pdPatchBoundary()")
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

//     Pout << "patch name: " << patchName_ << ", patch size: " << patch.size() << endl;

//     Random& rndGen = cloud_.rndGen();

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    if (isA<cyclicPolyPatch>(patch))
    {
        FatalErrorIn("pdCyclicBoundary::pdCyclicBoundary()")
            << "Patch: " << patchName_ << " is a cyclic boundary. It should be a patch." 
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }
    
//     totalPatchSurfaceArea_ = patchSurfaceArea_;
// 
//     if(Pstream::parRun())
//     {
// 	reduce(totalPatchSurfaceArea_, sumOp<scalar>);
//            
//     }
//     Pout << "faces: " << faces_ << endl;

    //- set the targeted fields
//     densities_.setSize(cells_.size(), scalar(0.0));
//     velocities_.setSize(cells_.size(), vector::zero);
//     temperatures_.setSize(cells_.size(), scalar(0.0));
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<pdPatchBoundary> pdPatchBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
{
    word pdPatchBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting pdPatchBoundaryModel "
         << pdPatchBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdPatchBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdPatchBoundary::New(const dictionary&) : " << endl
            << "    unknown pdPatchBoundary type "
            << pdPatchBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pdPatchBoundary>
	(
		cstrIter()(t, mesh, cloud, dict)
	);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdPatchBoundary::~pdPatchBoundary()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- needs to be called after the referred interaction list is built in pdCloud
void pdPatchBoundary::setBoundaryFields()
{
}

void pdPatchBoundary::measurePropertiesBeforeControl(pdParcel& p)
{

    if(measurePropertiesAtWall_)
    {
        label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        label wppLocalFace = wpp.whichFace(p.face());
    
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);
    
        const pdParcel::constantProperties& constProps(cloud_.constProps(p.typeId()));
    
        scalar m = constProps.mass();
    
        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);
    
        scalar U_dot_nw = p.U() & nw;
    
        vector Ut = p.U() - U_dot_nw*nw;
    
        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);
    
        cloud_.stdFields().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;
        cloud_.stdFields().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;
        cloud_.stdFields().linearKEBF()[wppIndex][wppLocalFace] += 0.5*m*(p.U() & p.U())*invMagUnfA;
        cloud_.stdFields().rotationalEBF()[wppIndex][wppLocalFace] += p.ERot()*invMagUnfA;
        cloud_.stdFields().rotationalDofBF()[wppIndex][wppLocalFace] += constProps.rotationalDegreesOfFreedom()*invMagUnfA;
        cloud_.stdFields().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
    
        // pre-interaction energy
        preIE_ = 0.5*m*(p.U() & p.U()) + p.ERot() + p.EVib();
    
        // pre-interaction momentum
        preIMom_ = m*p.U();
    }
}

void pdPatchBoundary::measurePropertiesAfterControl(pdParcel& p)
{
    if(measurePropertiesAtWall_)
    {
        label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        label wppLocalFace = wpp.whichFace(p.face());
    
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);
    
        const scalar deltaT = mesh_.time().deltaTValue();
    
        const pdParcel::constantProperties& constProps(cloud_.constProps(p.typeId()));
    
        scalar m = constProps.mass();
    
        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);
    
        scalar U_dot_nw = p.U() & nw;
    
        vector Ut = p.U() - U_dot_nw*nw;
    
        scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, VSMALL);
    
        cloud_.stdFields().rhoNBF()[wppIndex][wppLocalFace] += invMagUnfA;
        cloud_.stdFields().rhoMBF()[wppIndex][wppLocalFace] += m*invMagUnfA;
        cloud_.stdFields().linearKEBF()[wppIndex][wppLocalFace] += 0.5*m*(p.U() & p.U())*invMagUnfA;
        cloud_.stdFields().rotationalEBF()[wppIndex][wppLocalFace] += p.ERot()*invMagUnfA;
        cloud_.stdFields().rotationalDofBF()[wppIndex][wppLocalFace] += constProps.rotationalDegreesOfFreedom()*invMagUnfA;
        cloud_.stdFields().momentumBF()[wppIndex][wppLocalFace] += m*Ut*invMagUnfA;
    
        // post-interaction energy
        scalar postIE = 0.5*m*(p.U() & p.U()) + p.ERot() + p.EVib();
    
        // post-interaction momentum
        vector postIMom = m*p.U();
    
        scalar deltaQ = cloud_.nParticle()*(preIE_ - postIE)/(deltaT*fA);
        vector deltaFD = cloud_.nParticle()*(preIMom_ - postIMom)/(deltaT*fA);
    
        cloud_.stdFields().qBF()[wppIndex][wppLocalFace] += deltaQ;
        cloud_.stdFields().fDBF()[wppIndex][wppLocalFace] += deltaFD;
    }
}

void pdPatchBoundary::writeTimeData
(
    const fileName& pathName,
    const word& fileName,
    const List< Pair<scalar> >& data
)
{
    OFstream timeFile(pathName/fileName+".raw");

    if (timeFile.good())
    {
        timeFile << data << endl;
    }
    else
    {
        FatalErrorIn("pdPatchBoundary::writeTimeData()")
            << "Cannot open file " << timeFile.name()
            << abort(FatalError);
    }
}

void pdPatchBoundary::writeTimeData
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

void pdPatchBoundary::writeTimeData
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
        FatalErrorIn("void pdPatchBoundary::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

void pdPatchBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("patchBoundaryProperties");
}



const labelList& pdPatchBoundary::controlPatch() const
{
    return faces_;
}

const labelList& pdPatchBoundary::controlZone() const
{
    return cells_;
}

const word& pdPatchBoundary::patchName() const
{
    return patchName_;
}

const label& pdPatchBoundary::patchId() const
{
    return patchId_;
}



const scalar& pdPatchBoundary::density() const
{
    return density_;
}

scalar& pdPatchBoundary::density()
{
    return density_;
}

const vector& pdPatchBoundary::velocity() const
{
    return velocity_;
}
vector& pdPatchBoundary::velocity()
{
    return velocity_;
}

const scalar& pdPatchBoundary::temperature() const
{
    return temperature_;
}

scalar& pdPatchBoundary::temperature()
{
    return temperature_;
}


const scalarField& pdPatchBoundary::densityField() const
{
    return densities_;
}

scalarField& pdPatchBoundary::densityField()
{
    return densities_;
}

const vectorField& pdPatchBoundary::velocityField() const
{
    return velocities_;
}
vectorField& pdPatchBoundary::velocityField()
{
    return velocities_;
}

const scalarField& pdPatchBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& pdPatchBoundary::temperatureField()
{
    return temperatures_;
}

const bool& pdPatchBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& pdPatchBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
