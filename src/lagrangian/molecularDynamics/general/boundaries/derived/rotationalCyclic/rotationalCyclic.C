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

#include "rotationalCyclic.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(rotationalCyclic, 0);

addToRunTimeSelectionTable(cyclicBoundary, rotationalCyclic, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
rotationalCyclic::rotationalCyclic
(
    Time& t,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cyclicBoundary(t, mesh, dict),
    propsDict_(dict.subDict(typeName + "Properties"))

{

    scalar tolerance = 0.1;

    if (propsDict_.found("tolerance"))
    {
        tolerance = readScalar(propsDict_.lookup("tolerance"));
    }

    if
    (
           (theta_ < (mathematicalConstant::pi - SMALL))
        || (theta_ > (mathematicalConstant::pi + SMALL))
    )
    {
        rotate_ = true;

        if (propsDict_.found("theta"))
        {
            theta_ = readScalar(propsDict_.lookup("theta"));
            theta_ *= (mathematicalConstant::pi/180.0);
        }

        rotationStartPoint_ = propsDict_.lookup("rotationAxisStartPoint");
        rotationEndPoint_ = propsDict_.lookup("rotationAxisEndPoint");


        rotationAxis_ = (rotationEndPoint_ - rotationStartPoint_)
                        /mag(rotationEndPoint_ - rotationStartPoint_);
        // mid-point
        rotationPt_ = rotationAxis_*0.5*mag(rotationEndPoint_ - rotationStartPoint_)
                     + rotationStartPoint_;

        Info << "rotation mid-point : " << rotationPt_ << endl;

//         rotationAxis_ = nA_ ^ nB_;

        vector rotAxisA = rotationAxis_;
        vector rotAxisB = -rotationAxis_;

        Info << "rotational axis A = " << rotAxisA
                << ", rotational axis B = " << rotAxisB
                << endl;

        Info << "angle between cyclic boundaries is = " << theta_
                << ", ... reading in rotationPoint."
                << endl;

//         rotationPt_ = propsDict_.lookup("rotationPoint");



        // switching off translation for now
        cyclicTranslationVector_ = vector::zero;

        tensor aCrossA = tensor
                        (
                            0.0, -rotAxisA.z(), rotAxisA.y(),
                            rotAxisA.z(), 0.0, -rotAxisA.x(),
                            -rotAxisA.y(), rotAxisA.x(), 0.0
                        );

        tensor RAB = I*cos(theta_) + sin(theta_)*aCrossA
             + (1.0-cos(theta_))*rotAxisA*rotAxisA;

        tensor aCrossB = tensor
                        (
                            0.0, -rotAxisB.z(), rotAxisB.y(),
                            rotAxisB.z(), 0.0, -rotAxisB.x(),
                            -rotAxisB.y(), rotAxisB.x(), 0.0
                        );

        tensor RBA = I*cos(theta_) + sin(theta_)*aCrossB
             + (1.0-cos(theta_))*rotAxisB*rotAxisB;

        //A -> B
        {
            vector pADash = (RBA & (pA_ - rotationPt_)) + rotationPt_;
    
            Info<< "midpoint on boundary A, " << pA_ 
                << ", rotated (trial attempt 1) = " << pADash 
                << ", comparing with midpoint on boundary B: " << pB_
                << endl;
    
            if(mag(pADash - pB_) < tolerance )
            {
                RAB_ = RBA;
                Info << " accepted." << endl;
            }
            else
            {
                pADash = (RAB & (pA_ - rotationPt_)) + rotationPt_;
    
                Info<< "midpoint on boundary A, " << pA_ 
                    << ", rotated (trial attempt 2) = " << pADash 
                    << ", comparing with midpoint on boundary B: " << pB_
                    << endl;
    
                if(mag(pADash - pB_) < tolerance )
                {
                    RAB_ = RAB;
                    Info << " accepted." << endl;
                }
                else
                {
                    FatalErrorIn("rotationalCyclic::rotationalCyclic()")
                        << "Rotation of boundary midpoint A: " << pA_ 
                        << " does not coindcide with midpoint on boundary B: " << pB_ 
                        << ". Please pick midpoints on both boundaries so that they coincide exactly. "
                        << " Change tolerance if necessary: " << tolerance 
                        << nl << "in: "
                        << t.system()/"boundariesDict"
                        << exit(FatalError);
                }
            }
        }

        //B -> A
        {
            vector pBDash = (RBA & (pB_ - rotationPt_)) + rotationPt_;
    
            Info<< "midpoint on boundary B, " << pB_ 
                << ", rotated (trial attempt 1) = " << pBDash 
                << ", comparing with midpoint on boundary A: " << pA_
                << endl;
    
            if(mag(pBDash - pA_) < tolerance )
            {
                RBA_ = RBA;
                Info << " accepted." << endl;
            }
            else
            {
                pBDash = (RAB & (pB_ - rotationPt_)) + rotationPt_;
    
                Info<< "midpoint on boundary B, " << pB_ 
                    << ", rotated (trial attempt 2) = " << pBDash 
                    << ", comparing with midpoint on boundary A: " << pA_
                    << endl;
    
                if(mag(pBDash - pA_) < tolerance )
                {
                    RBA_ = RAB;
                    Info << " accepted." << endl;
                }
                else
                {
                    FatalErrorIn("rotationalCyclic::rotationalCyclic()")
                        << "Rotation of boundary midpoint B: " << pB_ 
                        << " does not coindcide with midpoint on boundary A: " << pA_ 
                        << ". Please pick midpoints on both boundaries so that they coincide exactly. "
                        << " Change tolerance if necessary: " << tolerance 
                        << nl << "in: "
                        << t.system()/"boundariesDict"
                        << exit(FatalError);
                }
            }
        }
    }
    else
    {
        FatalErrorIn("rotationalCyclic::rotationalCyclic()")
            << "Patch: " << patchName_ << " does not need a rotational cyclic boundary model. " 
            << " Angle = " << theta_
            << nl << "in: "
            << t.system()/"boundariesDict"
            << exit(FatalError);
    }

//     Info << "TESTING" << endl;

    //serial - not parallelised
//     labelList coupledFacesA = coupledFacesA_;
//     labelList coupledFacesB = coupledFacesB_;
//     labelList cellsB = cellsB_;
// 
// //     Info << "coupledFacesA: "  << coupledFacesA << endl;
// //     Info << "coupledFacesB: "  << coupledFacesB << endl;
//     label nCoupled = 0;
//     forAll(coupledFacesA, fA)
//     {
//         const label& faceA = coupledFacesA[fA];
//         const vector& fCA = mesh_.faceCentres()[faceA];
// 
//         forAll(coupledFacesB, fB)  
//         {
//             const label& faceB = coupledFacesB[fB];
//             const vector& fCB = mesh_.faceCentres()[faceB];
// 
//             vector fCBDash = (RBA_ & (fCB - rotationPt_)) + rotationPt_;
// //             Info << "face centre: " << fCB << " trial face matching centre: "
// //                  << fCBDash << endl;
//             if(mag(fCBDash - fCA) < tolerance)
//             {
//                 coupledFacesB_[fA] = faceB;
//                 cellsB_[fA] = cellsB[fB];
//                 nCoupled++;
//             }
//         }
//     }

//     Info << "nCoupled: "<< nCoupled << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rotationalCyclic::~rotationalCyclic()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
