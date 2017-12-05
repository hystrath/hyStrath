/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    polyFaceTracker

Description

\*----------------------------------------------------------------------------*/

#include "polyFaceTracker.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "polyMoleculeCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
polyFaceTracker::polyFaceTracker
(
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    molCloud_(molCloud),
    molIdFlux_(molCloud_.cP().molIds().size()),
    massIdFlux_(molCloud_.cP().molIds().size()),
    absMomIdFlux_(molCloud_.cP().molIds().size()),
    momIdFlux_(molCloud_.cP().molIds().size())

{
    forAll(molIdFlux_, i)
    {
        molIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
        massIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
        absMomIdFlux_[i].setSize(mesh_.nFaces(), 0.0);
        momIdFlux_[i].setSize(mesh_.nFaces(), vector::zero);
    }

//     Info << "mesh boundary names = " << mesh_.boundaryMesh().names() << endl;

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyFaceTracker::~polyFaceTracker()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyFaceTracker::clean()
{
    //- clean geometric fields

    forAll(molIdFlux_, i)
    {
        molIdFlux_[i] = scalar(0.0);
        massIdFlux_[i] = scalar(0.0);
        absMomIdFlux_[i] = scalar(0.0);
        momIdFlux_[i] = vector::zero;
    }
}



// called during the move function 
// 
void polyFaceTracker::updateFields
(
    polyMolecule& mol
)
{
    const label& crossedFace = mol.face();
    const label& molId = mol.id();
//     const polyMolecule::constantProperties& constProp = molCloud_.constProps(molId);
    
    const scalar& mass = molCloud_.cP().mass(molId);
    const vector& U = mol.v();
    const vector mom = mol.v()*mass;
//     const scalar pE = mol.potentialEnergy();
//     const scalar kE = 0.5*mass*magSqr(U);
//     const scalar energy = kE + pE;
//     const vector force = mol.a()*mass;
//     const scalar volume = 4.0*mathematicalConstant::pi*mol.R()*mol.R()*mol.R()/3.0;


    const label& patchId = mesh_.boundaryMesh().whichPatch(crossedFace);

    vector nF = mesh_.faceAreas()[crossedFace];
    nF /= mag(nF);

//     Pout << "Mol at pos = " << mol.position() 
//          << ", velocity = " << mol.v()
//          << ", face normal = " << nF
//          << ", face centre = " <<  mesh_.faceCentres()[crossedFace]
//          << endl;
         
    //- direction of polyMolecule trajectory with respect to the face normal
    scalar sgn = sign( U & mesh_.faceAreas()[crossedFace] ) * 1.0;

    //- geometric fields

    if(patchId != -1) //- boundary face
    {
/*        Pout<< "patchId = " << patchId
            << ", patchName = " << mesh_.boundaryMesh().names()[patchId]
            << endl;  */      
        
        const polyPatch& patch = mesh_.boundaryMesh()[patchId];

        const label faceIndex = crossedFace - patch.start();

        //- cyclic patches
        // NOTE: in a purely cyclic condition, a molecule is passed to the coupled *receiving* face 
        //       before calling this function.
        if (isA<cyclicPolyPatch>(patch))
        {
//             Pout << "Cyclic" << endl;
            
            label coupledFace = refCast<const cyclicPolyPatch>
            (
                patch
            ).neighbPatch().start() + faceIndex;

            molIdFlux_[molId][coupledFace] += 1.0;
            massIdFlux_[molId][coupledFace] += mass;
            absMomIdFlux_[molId][crossedFace] += mag(mom & nF);
            momIdFlux_[molId][crossedFace] += mom;            
        }
        
        //- processor patches
        // NOTE: properties are appended to the face of the *leaving* processor only.
        // normal vector points out from the domain.
        // NOTE: when the processor patch is also a cyclic patch,
        // properties still remain the same -> i.e. *leaving* processor face
        if (isA<processorPolyPatch>(patch))
        {
//             Pout << "Processor" << endl;
            
            molIdFlux_[molId][crossedFace] += sgn*1.0;
            massIdFlux_[molId][crossedFace] += sgn*mass;
            absMomIdFlux_[molId][crossedFace] += mag(mom & nF);
            momIdFlux_[molId][crossedFace] += mom;
        }
    }
    else //- internal face
    {
//         Info << "internal" << endl;
        
        molIdFlux_[molId][crossedFace] += sgn*1.0;
        massIdFlux_[molId][crossedFace] += sgn*mass;
        absMomIdFlux_[molId][crossedFace] += mag(mom & nF);
        momIdFlux_[molId][crossedFace] += mom;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
