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
    constructCyclicBoundaries

Description

\*----------------------------------------------------------------------------*/

#include "constructCyclicBoundaries.H"
#include "polyBoundaryMeshEntries.H"
#include "graph.H"
#include "cyclicPolyPatch.H"
namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


constructCyclicBoundaries::constructCyclicBoundaries
(
    const polyMesh& mesh
)
:
    mesh_(mesh)
{
    Info << "Constructing information for cyclic boundaries"
         << endl;

    const polyBoundaryMesh& bM = mesh_.boundaryMesh();

    nCyclicBoundaries_ = 0;

    forAll(bM, patchI)
    {
        const polyPatch& patch = bM[patchI];
        
        if (isA<cyclicPolyPatch>(patch))
        {
            nCyclicBoundaries_++;
        }
    }
    
    Info << "Number of cyclic boundaries = " << nCyclicBoundaries_<< endl;

    boundaryPoints_.setSize(nCyclicBoundaries_);
    names_.setSize(nCyclicBoundaries_);    
    nFaces_.setSize(nCyclicBoundaries_);
    normal_.setSize(nCyclicBoundaries_);
    centroid_.setSize(nCyclicBoundaries_);
    neighbourNames_.setSize(nCyclicBoundaries_);
    label c = 0;

    forAll(bM, patchI)
    {
        const polyPatch& patch = bM[patchI];
        
        if (isA<cyclicPolyPatch>(patch))
        {
            names_[c] = bM.names()[patchI];
            nFaces_[c] = patch.size();
            const cyclicPolyPatch& cyclicPatch =
                            refCast<const cyclicPolyPatch>(patch);
                            
            neighbourNames_[c] = cyclicPatch.neighbPatchName();            
            setPatchPoints(patchI, c);
            
            c++;
        }
    }

    setCyclicBoundaries();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constructCyclicBoundaries::~constructCyclicBoundaries()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




void constructCyclicBoundaries::setPatchPoints
(
    const label& patchId,
    const label& c
)
{
    const polyBoundaryMesh& bM = mesh_.boundaryMesh();
    const polyPatch& patch = bM[patchId];

    DynamicList<label> labelPoints;
    DynamicList<vector> points;

    label nFaces = 0;
    scalar A = 0.0;
    vector nF = vector::zero;    
    vector CA = vector::zero;
    
    for(label i = 0; i < patch.size(); i++)
    {
        nFaces++;            
        
        label globalFaceI = patch.start() + i;
        
        // normal surface area vector
        vector sF = mesh_.faceAreas()[globalFaceI]; 
        
        // area of face
        scalar area = Foam::mag(sF);
        
        A += area; 
        nF += sF/area; // normal vector
        const vector& C = mesh_.faceCentres()[globalFaceI];// centroid
        CA += C*area;
        
        forAll(mesh_.faces()[globalFaceI], j)
        {
            const label& p = mesh_.faces()[globalFaceI][j];
            
            if(findIndex(labelPoints, p) == -1)
            {
                labelPoints.append(p);
                points.append(mesh_.points()[p]);
            }
        }
    }
    
//     Info << nl << "boundary = " << names_[c] << endl;    
    
   
    //points.shrink();

    normal_[c] = nF/nFaces;
    
    centroid_[c] = CA/A;

//     Info<< "centroid = " << centroid_[c]
//         << ", normal = " << normal_[c]
//         << endl;
        
    // find the bounding 4 points
    
    vector pointA = vector::zero;
    vector pointB = vector::zero;
    vector pointC = vector::zero;    
    vector pointD = vector::zero;

    // a - find the point which is farthest away from centroid 
    //   - if a rectangular plane, it is likely to be one of the corner points
    
    scalar rMag = 0.0;
    
    forAll(points, i)
    {
        scalar rD = mag(points[i]-centroid_[c]);
        
        if(rD > rMag)
        {
            rMag = rD;
            pointA = points[i];
        }
    }
    
//     Info << "first point = " << pointA << endl;
    
    // b- second point can be found along the diagonal of the rectangular point.
    
    vector refV = -(pointA-centroid_[c])/rMag;
    pointB = centroid_[c] + rMag*refV;
    
//     Info << "second point (est) = " << pointB << endl;    

    //- check
    
    scalar tolerance = 0.1;
    
    forAll(points, i)
    {
        scalar rD = mag(points[i]-pointB);
        
//         Info << "point = " << points[i] << ", rD = " << rD << endl;
        
        if(rD < tolerance)
        {
            pointB = vector::zero;
            pointB = points[i];
//             Info << "second point found = " << points[i] << endl;
        }
    }
    
//     Info << "second point = " << pointB << endl;    
    
    // c, - third point is found by now starting from point B
    //   and finding the farthest point in a third direction.
    
    vector refV2 = normal_[c] ^ (-refV);
    scalar rMag2 = 0.0;
    
    forAll(points, i)
    {
        scalar rD = mag((points[i]-pointB) & refV2);
        
        if(rD > rMag2)
        {
            rMag2 = rD;
            pointC = points[i];
        }
    }    
    
//     Info << "third point = " << pointC << endl;    
    
    // d. fourth point is found by repeating step b.
    vector refV3 = -(pointC-centroid_[c])/rMag2;
    pointD = centroid_[c] + rMag2*refV3;
    
//     Info << "fourth point (est) = " << pointD << endl;    
    // check
    
    forAll(points, i)
    {
        scalar rD = mag(points[i]-pointD);
        
        if(rD < tolerance)
        {
            pointD = vector::zero;            
            pointD = points[i];
//             Info << "fourth point found = " << points[i] << endl;            
        }
    }   
    
//     Info << "fourth point = " << pointD << endl;      
    
    boundaryPoints_[c].setSize(4, vector::zero);
    
    boundaryPoints_[c][0]=pointA;
    boundaryPoints_[c][1]=pointD;
    boundaryPoints_[c][2]=pointB;
    boundaryPoints_[c][3]=pointC;

    
    // INFO
    
    Info << nl << "Cyclic Boundary Name = " << names_[c]
         << ", Neighbour patch Name = " << neighbourNames_[c]
         << endl;       
    
    Info<< nl << "centroid = " << centroid_[c]
        << nl << "normal = " << normal_[c]
        << nl << "nFaces = " << nFaces_[c]
        << nl << "bounding Pts = " << boundaryPoints_[c]
        << endl;
}


void constructCyclicBoundaries::setCyclicBoundaries()
{
    Info << nl <<  "Output" << nl << endl;

    // create file names 
    fileName pathName = mesh_.time().path()/"system";
    fileName nameFile = "cyclicBoundaries";

    {
        // deletes current content of file
        OFstream file(pathName/nameFile);
    
        if(file.good())
        {
            file << endl;
        }
        else
        {
            FatalErrorIn("void writeTimeData::writeTimeData()")
                << "Cannot open file " << file.name()
                << abort(FatalError);
        }
    }

    // creates a new file
    IOdictionary dict
    (
        IOobject
        (
            "cyclicBoundaries",
            mesh_.time().system(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    dict.clear();

    if (!mesh_.write())     // - requrired for zone writing
    {
            FatalErrorIn("constructCyclicBoundaries::setCyclicBoundaries()")
                << "Failed writing mesh."
                << exit(FatalError);
    }

    // appends to the file
    {
        fileName fName(pathName/nameFile);
    
        std::ofstream file(fName.c_str(),ios_base::app);
    
        if(file.is_open())
        {
            file << nl << "cyclicBoundaries"  << nl
                    << "(" << nl;

    
            const polyBoundaryMesh& bM = mesh_.boundaryMesh();

            label c = 0;

            forAll(bM, patchI)
            {
                const polyPatch& patch = bM[patchI];
        
                if (isA<cyclicPolyPatch>(patch))
                {
                    file << "\t"<< "boundary"  << nl
                         << "\t" << "{" << nl
                         << "\t" << "\t" << "cyclicBoundaryProperties" << nl
                         << "\t" << "\t" << "{" << nl
                         << "\t" << "\t" << "\t" << "patchName" 
                         << "\t" << "\t" << "\t" << names_[c] << ";" << nl
                         << "\t" << "\t" << "\t" << "neighbourPatchName" 
                         << "\t" << "\t" << "\t" << neighbourNames_[c] << ";" << nl
                         << "\t" << "\t" << "\t" << "nFaces" 
                         << "\t" << "\t" << "\t" << nFaces_[c] << ";" << nl
                         << "\t" << "\t" << "\t" << "normal" 
                         << "\t" << "\t" << "\t" << "(" << normal_[c].x() 
                         << " "<< normal_[c].y() << " " << normal_[c].z() << ");" << nl
                         << "\t" << "\t" << "\t" << "boundaryPoints" << nl
                         << "\t" << "\t" << "\t" << "(" << nl;

                    forAll(boundaryPoints_[c], p)
                    {
                        file << "\t" << "\t" << "\t" << "\t"<< "(" << boundaryPoints_[c][p].x() 
                             << " "<< boundaryPoints_[c][p].y() << " " << boundaryPoints_[c][p].z()
                             << ")" << nl;
                    }

                    file << "\t" << "\t" << "\t" << ");" << nl
                         << "\t" << "\t" << "}" << nl
                         << "\t" << "\t" <<  "cyclicBoundaryModel" 
                         << "\t" << "\t" <<  "standardCyclic" << ";" << nl
                         << "\t"  << "}" << nl;

                    c++;
                }
            }
                    
            file << ");" << nl;

            file<< nl 
                << "// ************************************************************************* //" 
                << nl;
        }
        else
        {
            FatalErrorIn("void constructCyclicBoundaries::constructCyclicBoundaries()")
                << "Cannot open file " << fName
                << abort(FatalError);
        }
    
        file.close();
    }
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
