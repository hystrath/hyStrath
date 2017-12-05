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

#include "polyStandardCyclic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyStandardCyclic, 0);

addToRunTimeSelectionTable(polyCyclicBoundary, polyStandardCyclic, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyStandardCyclic::polyStandardCyclic
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyCyclicBoundary(t, mesh, molCloud, dict)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyStandardCyclic::~polyStandardCyclic()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyStandardCyclic::calculateProperties()
{}

void polyStandardCyclic::initialConfiguration()
{}

void polyStandardCyclic::controlAfterMove()
{}

void polyStandardCyclic::controlMol
(
    polyMolecule& mol,
    polyMolecule::trackingData& td,
    const label& patchi
)
{
    /*
    label& faceI = mol.face();
//     Info << "Position: " << mol.position() << endl;
//     Info<< "Before hit. faceI " << faceI 
//         << ", fc: " << mesh_.faceCentres()[faceI]
//         << ", cellI: " << mol.cell()
//         << endl;

//     Pout << "error1: faceI: " << faceI << endl;
//     Pout << " error1: fC: " << mesh_.faceCentres()[faceI]<< endl;

//     label patchI = mesh_.boundaryMesh().whichPatch(faceI);
//     Pout << "error1: patchi: " << patchi << endl;
    
    // hit cyclic patch (serial)
    if(cyclic()->patchId() == patchi)
    {
//         Pout <<"hit cyclic patch: " << mol.position() << endl;
//         Info << "hit poly patch name: " << patchName_ << endl;
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        const cyclicPolyPatch& cpp = static_cast<const cyclicPolyPatch&>(patch);

        // modifying faceI and cellI
        faceI = cpp.transformGlobalFace(faceI);

        mol.cell() = mesh_.faceOwner()[faceI];

//         Info<< "modify face " << faceI 
//             << ", fc: " <<  mesh_.faceCentres()[faceI]
//             << ", cell: " << mol.cell()
//             << endl;

        // move particle
        const vector& nF = mesh_.faceAreas()[faceI];
        mol.position() += nF/mag(nF)*mag(cyclic()->cyclicTranslationVector());

//         Pout << "Position after hit: " << mol.position() << endl;
    }
    else if(patchi != -1)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        // hit processor patch (in parallel) but it is a cyclic patch
        if( isA<processorPolyPatch>(patch))
        {
//             Pout << "error3" << endl;
            //belongs to this cylic 
            if(findIndex(cyclic()->controlPatch(), faceI) != -1)
            {
//                 Pout <<"hit processor patch: " << mol.position() << endl;
    
                const vector& nF = mesh_.faceAreas()[faceI];
                mol.position() -= nF/mag(nF)*mag(cyclic()->cyclicTranslationVector());
    
//                 Pout <<"new position: " << mol.position() << endl;
            }
        }
    }
//     Pout << "error2" << endl;

//     Pout << "hit cyclic boundary (new model)! "<< mol.position() <<  endl;
    


//     Pout << "hit cyclic boundary (new model)! "<< mol.position() <<  endl;

//     if(!isA<processorPolyPatch>(patch))
//     {
// //         Pout << "hit cyclic boundary (new model)! "<< mol.position() <<  endl;
//     
//         const vector& nF = mesh_.faceAreas()[faceI];
//     
//         mol.position() -= nF/mag(nF)*mag(cyclicTranslationVector_);
//     
// //         Pout << "hit cyclic boundary (after new model)! "<< mol.position() <<  endl;
//     }

//     const vector& fC  = mesh_.faceCentres()[faceI];

//     label patchFacei_ = cpp.whichFace(facei_);

//     const label f = findIndex(controlPatch(), faceI);

//     label patchi = mesh_.boundaryMesh().whichPatch(faceI);

//     Info << "mol at pos: " << mol.position() 
//         << ", new pos: " << mol.position() - nF/mag(nF)*mag(cyclicTranslationVector_) 
//          << ", v = " << mol.v() 
//          << " hits patch name: " << mesh_.boundaryMesh().names()[patchi]
//          << endl;
    
   
//     const polyPatch& patch = mesh_.boundaryMesh()[patchi];

//     const cyclicPolyPatch& cpp = static_cast<const cyclicPolyPatch&>(patch);

//     Info << "patch is not parallel: " << !cpp.parallel() << endl;
//     Info << "patch is separated: " << cpp.separated() << endl;


//     label patchFacei_ = cpp.whichFace(facei_);

//     facei_ = cpp.transformGlobalFace(facei_);

//     celli_ = cloud_.polyMesh_.faceOwner()[facei_];

//     if (!cpp.parallel())
//     {
//         const tensor& T = cpp.transformT(patchFacei_);
// 
//         transformPosition(T);
//         static_cast<ParticleType&>(*this).transformProperties(T);
//     }
//     else if (cpp.separated())
//     {
//         position_ += cpp.separation(patchFacei_);
//         static_cast<ParticleType&>(*this).transformProperties
//         (
//             cpp.separation(patchFacei_)
//         );
//     }

// 
// 
//     label patchFacei_ = cpp.whichFace(facei_);
// 
//     facei_ = cpp.transformGlobalFace(facei_);
// 
//     celli_ = cloud_.polyMesh_.faceOwner()[facei_];
// 
//     if (!cpp.parallel())
//     {
//         const tensor& T = cpp.transformT(patchFacei_);
// 
//         transformPosition(T);
//         static_cast<ParticleType&>(*this).transformProperties(T);
//     }
//     else if (cpp.separated())
//     {
//         position_ += cpp.separation(patchFacei_);
//         static_cast<ParticleType&>(*this).transformProperties
//         (
//             cpp.separation(patchFacei_)
//         );
//     }
*/
}

// void polyStandardCyclic::correctAfterParallelTransfer
// (
//     polyMolecule& mol,
//     const label& patchi,
//     polyMolecule::trackData& td
// )
// {
// 
// }


void polyStandardCyclic::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}


void polyStandardCyclic::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

}




} // End namespace Foam

// ************************************************************************* //
