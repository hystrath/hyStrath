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

#include "dsmcPatchBoundary.H"
#include "IFstream.H"
#include "graph.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcPatchBoundary, 0);

defineRunTimeSelectionTable(dsmcPatchBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcPatchBoundary::dsmcPatchBoundary
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
    boundaryDict_(dict.subDict("patchBoundaryProperties")),
    patchName_(boundaryDict_.lookup("patchName")),
    patchId_(0),
    faces_(),
    nFaces_(0),
    patchSurfaceArea_(0.0),
    totalPatchSurfaceArea_(0.0),
    cells_(),
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
        FatalErrorIn("dsmcPatchBoundary::dsmcPatchBoundary()")
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    if (isA<cyclicPolyPatch>(patch))
    {
        FatalErrorIn("dsmcCyclicBoundary::dsmcCyclicBoundary()")
            << "Patch: " << patchName_ 
            << " is a cyclic boundary. It should be a patch." 
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
        
        //- area on one processor
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]); 
    }
    
    totalPatchSurfaceArea_ = patchSurfaceArea_;

    if(Pstream::parRun())
    {
        //- total area on all processors
        reduce(totalPatchSurfaceArea_, sumOp<scalar>());  
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcPatchBoundary> dsmcPatchBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcPatchBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting dsmcPatchBoundaryModel "
         << dsmcPatchBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcPatchBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcPatchBoundary::New(const dictionary&) : " << endl
            << "    unknown dsmcPatchBoundary type "
            << dsmcPatchBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid patch boundary types are :" << endl;
            
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcPatchBoundary>
    (
        cstrIter()(t, mesh, cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcPatchBoundary::~dsmcPatchBoundary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- need to be called after the referred interaction list is built in dsmcCloud
void dsmcPatchBoundary::setBoundaryFields()
{}


void dsmcPatchBoundary::setNewBoundaryFields()
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    
    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());
    
    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    //- loop through all faces and set the boundary cells
    //  no conflict with parallelisation because the faces are unique
    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        
        //- area on one processor
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]); 
    }
    
    totalPatchSurfaceArea_ = patchSurfaceArea_;

    if(Pstream::parRun())
    {
        //- total area on all processors
        reduce(totalPatchSurfaceArea_, sumOp<scalar>());  
    }
}


void dsmcPatchBoundary::calculateWallUnitVectors
(
    dsmcParcel& p,
    vector& nw, 
    vector& tw1, 
    vector& tw2
)
{
    nw = p.normal();
    nw /= mag(nw);
    
    vector& U = p.U();

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;
    
    Random& rndGen = cloud_.rndGen();
    
    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    tw2 = nw^tw1;
}


void dsmcPatchBoundary::measurePropertiesBeforeControl(dsmcParcel& p)
{
    if(measurePropertiesAtWall_)
    {
        const label& wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());

        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);
        
        const scalar deltaT = cloud_.deltaTValue(p.cell());
    
        const dsmcParcel::constantProperties& 
            constProps(cloud_.constProps(p.typeId()));

        const scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        const scalar U_dot_nw = p.U() & nw;
    
        const vector& Ut = p.U() - U_dot_nw*nw;
    
        const scalar invMagUnfADt = 1.0/max(mag(U_dot_nw)*fA*deltaT, SMALL);

        //- Update boundary flux measurements
        cloud_.boundaryFluxMeasurements()
            .rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt;
        
        if(constProps.rotationalDegreesOfFreedom() > 0)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNIntBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt;
        }

        if(constProps.numberOfElectronicLevels() > 1)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNElecBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt; 
        }

        cloud_.boundaryFluxMeasurements()
            .rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*invMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .linearKEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                0.5*m*(p.U() & p.U())*invMagUnfADt;
            
        cloud_.boundaryFluxMeasurements()
            .mccSpeciesBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*(p.U() & p.U())*invMagUnfADt;
            
        cloud_.boundaryFluxMeasurements()
            .momentumBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*Ut*invMagUnfADt;
            
        cloud_.boundaryFluxMeasurements()
            .rotationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                p.ERot()*invMagUnfADt;
            
        cloud_.boundaryFluxMeasurements()
            .rotationalDofBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.rotationalDegreesOfFreedom()*invMagUnfADt;

        forAll(p.vibLevel(), mode)
        {
            cloud_.boundaryFluxMeasurements()
                .vibrationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                    p.vibLevel()[mode]*constProps.thetaV()[mode]
                  * physicoChemical::k.value()
                  * invMagUnfADt;
                  
            cloud_.boundaryFluxMeasurements()
                .evmsBF()[p.typeId()][mode][wppIndex][wppLocalFace] += 
                    p.vibLevel()[mode]*constProps.thetaV()[mode]
                  * physicoChemical::k.value()
                  * invMagUnfADt;      
        }

        cloud_.boundaryFluxMeasurements()
            .electronicEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.electronicEnergyList()[p.ELevel()]
              * invMagUnfADt;

        //- pre-interaction energy
        preIE_ = 0.5*m*(p.U() & p.U()) + p.ERot() 
            + constProps.electronicEnergyList()[p.ELevel()];

        forAll(p.vibLevel(), mode)
        {
           preIE_ += p.vibLevel()[mode]*constProps.thetaV()[mode]
              *physicoChemical::k.value();
        }

        //- pre-interaction momentum
        preIMom_ = m*p.U();
    }
}

void dsmcPatchBoundary::measurePropertiesAfterControl
(
    dsmcParcel& p, 
    const scalar& heatOfReaction
)
{
    if(measurePropertiesAtWall_)
    {       
        const label wppIndex = patchId_;
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = wpp.whichFace(p.face());
    
        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);
    
        const scalar deltaT = cloud_.deltaTValue(p.cell());
    
        const dsmcParcel::constantProperties& 
            constProps(cloud_.constProps(p.typeId()));
    
        const scalar m = constProps.mass();
    
        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);
    
        const scalar U_dot_nw = p.U() & nw;
    
        const vector Ut = p.U() - U_dot_nw*nw;
    
        const scalar invMagUnfADt = 1.0/max(mag(U_dot_nw)*fA*deltaT, SMALL);

        //- Update boundary flux measurements
        cloud_.boundaryFluxMeasurements()
            .rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt;
        
        if(constProps.rotationalDegreesOfFreedom() > 0)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNIntBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt; 
        }
        
        if(constProps.numberOfElectronicLevels() > 1)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNElecBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfADt;
        }
        
        cloud_.boundaryFluxMeasurements()
            .rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*invMagUnfADt;
            
        cloud_.boundaryFluxMeasurements()
            .linearKEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                0.5*m*(p.U() & p.U())*invMagUnfADt;
        
        cloud_.boundaryFluxMeasurements()
            .mccSpeciesBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*(p.U() & p.U())*invMagUnfADt;
        
        cloud_.boundaryFluxMeasurements()
            .momentumBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*Ut*invMagUnfADt;
        
        cloud_.boundaryFluxMeasurements()
            .rotationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                p.ERot()*invMagUnfADt;
        
        cloud_.boundaryFluxMeasurements()
            .rotationalDofBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.rotationalDegreesOfFreedom()*invMagUnfADt;
        
        forAll(p.vibLevel(), mode)
        {
            cloud_.boundaryFluxMeasurements()
                .vibrationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                    p.vibLevel()[mode]*constProps.thetaV()[mode]
                  * physicoChemical::k.value()
                  * invMagUnfADt;
                  
            cloud_.boundaryFluxMeasurements()
                .evmsBF()[p.typeId()][mode][wppIndex][wppLocalFace] += 
                    p.vibLevel()[mode]*constProps.thetaV()[mode]
                  * physicoChemical::k.value()
                  * invMagUnfADt;       
        }
        
        cloud_.boundaryFluxMeasurements()
            .electronicEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.electronicEnergyList()[p.ELevel()]*invMagUnfADt;
        
        //- post-interaction energy
        scalar postIE = 0.5*m*(p.U() & p.U()) + p.ERot() 
            + constProps.electronicEnergyList()[p.ELevel()];
        
        forAll(p.vibLevel(), mode)
        {
           postIE += p.vibLevel()[mode]*constProps.thetaV()[mode]
              *physicoChemical::k.value();
        }
        
        //- post-interaction momentum
        const vector& postIMom = m*p.U();
    
        const scalar nParticle = cloud_.nParticles(wppIndex, wppLocalFace);
        
        const scalar deltaQ = nParticle
            * (preIE_ - postIE + (heatOfReaction*physicoChemical::k.value()))
            / (deltaT*fA);
            
        const vector deltaFD = nParticle*(preIMom_ - postIMom)/(deltaT*fA);
        
        cloud_.boundaryFluxMeasurements()
            .qBF()[p.typeId()][wppIndex][wppLocalFace] += deltaQ;
            
        cloud_.boundaryFluxMeasurements()
            .fDBF()[p.typeId()][wppIndex][wppLocalFace] += deltaFD;
    }
}


void dsmcPatchBoundary::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData
)
{
    fileName writeFile(pathName/nameFile);

    graph outputGraph("title", "x", "y", xData);

    outputGraph.write(writeFile, "raw");
}


void dsmcPatchBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("patchBoundaryProperties");
}


const labelList& dsmcPatchBoundary::controlPatch() const
{
    return faces_;
}


const labelList& dsmcPatchBoundary::controlZone() const
{
    return cells_;
}


const word& dsmcPatchBoundary::patchName() const
{
    return patchName_;
}


const label& dsmcPatchBoundary::patchId() const
{
    return patchId_;
}


const scalar& dsmcPatchBoundary::density() const
{
    return density_;
}


scalar& dsmcPatchBoundary::density()
{
    return density_;
}


const vector& dsmcPatchBoundary::velocity() const
{
    return velocity_;
}


vector& dsmcPatchBoundary::velocity()
{
    return velocity_;
}


const scalar& dsmcPatchBoundary::temperature() const
{
    return temperature_;
}


scalar& dsmcPatchBoundary::temperature()
{
    return temperature_;
}


const scalarField& dsmcPatchBoundary::densityField() const
{
    return densities_;
}


scalarField& dsmcPatchBoundary::densityField()
{
    return densities_;
}


const vectorField& dsmcPatchBoundary::velocityField() const
{
    return velocities_;
}


vectorField& dsmcPatchBoundary::velocityField()
{
    return velocities_;
}


const scalarField& dsmcPatchBoundary::temperatureField() const
{
    return temperatures_;
}


scalarField& dsmcPatchBoundary::temperatureField()
{
    return temperatures_;
}


const bool& dsmcPatchBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}


const bool& dsmcPatchBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
