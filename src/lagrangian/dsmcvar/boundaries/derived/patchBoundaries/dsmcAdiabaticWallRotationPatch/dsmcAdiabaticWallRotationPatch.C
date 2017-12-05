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

#include "dsmcAdiabaticWallRotationPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAdiabaticWallRotationPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcAdiabaticWallRotationPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAdiabaticWallRotationPatch::dsmcAdiabaticWallRotationPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    wallVelocity_ = readScalar(propsDict_.lookup("velocity"));
    referenceTemperature_ = readScalar(propsDict_.lookup("referenceTemperature"));
    rotationAxis_ = propsDict_.lookup("rotationAxis");
    centrePoint_ = propsDict_.lookup("centrePoint");
    rotationAxis_ /= mag(rotationAxis_);
    
    
    // test

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAdiabaticWallRotationPatch::~dsmcAdiabaticWallRotationPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcAdiabaticWallRotationPatch::initialConfiguration()
{}

void dsmcAdiabaticWallRotationPatch::calculateProperties()
{

}

void dsmcAdiabaticWallRotationPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    
    // wall velocity
    vector uNew = (
                    (p.position() - centrePoint_)/mag((p.position() - centrePoint_))
                  )
                  ^ rotationAxis_;
                  
    uNew /= mag(uNew);
    uNew *= wallVelocity_;
    
    vector& U = p.U();
    
    measurePropertiesBeforeControl(p);
    
    scalar EInc = magSqr(U - uNew);
    
    U -= uNew;
    
//     Info << "Energy before = " << magSqr(U) << endl;
    
    
//     Info << "Energy before = " << EInc << endl;

    label typeId = p.typeId();

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    Random& rndGen(cloud_.rndGen());

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.scalar01()),
            U.y()*(0.8 + 0.2*rndGen.scalar01()),
            U.z()*(0.8 + 0.2*rndGen.scalar01())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

//         scalar T = boundaryT_.boundaryField()[wppIndex][wppLocalFace];

    scalar mass = cloud_.constProps(typeId).mass();   
    
//     Info << "Position = " << p.position() << endl;
//     Info << "Velocity = " << uNew << endl;

//     while(true)
//     {
//         U =
//             sqrt(physicoChemical::k.value()*T/mass)
//             *(
//                 rndGen.GaussNormal()*tw1
//                 + rndGen.GaussNormal()*tw2
//                 - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
//             );
//         
//         vector Urefl = U;
// 
//         scalar reScale = 1.0;
//         
//         //quadratic formula
//         
//         scalar a = magSqr(Urefl);
//         
//         scalar b = 2.0*(Urefl & uNew);
//             
//         scalar c = magSqr(uNew) - magSqr(Uinc);
//         
//         if( (sqr(b) - 4.0*a*c) > VSMALL)
//         {
//             reScale = (-b + sqrt(sqr(b) - 4.0*a*c))/(2.0*a);
//             
//             U *= reScale;
//                 
//             if( (((U & nw) + (nw & uNew)))/(Urefl & nw) < VSMALL)
//             {
//                 reScale = (-b - sqrt(sqr(b) - 4.0*a*c))/(2.0*a); 
//                 
//                 U = Urefl*reScale;
//             }
// 
//             if( (((U & nw) + (nw & uNew)))/(Urefl & nw) > VSMALL)
//             {   
//                 U += uNew;
// //                 Info << "Energy final = " << magSqr(U) << endl;
//                 break;
//             }
//         }
//     }
    
    U =
        sqrt(physicoChemical::k.value()*referenceTemperature_/mass)
        *(
            rndGen.GaussNormal()*tw1
            + rndGen.GaussNormal()*tw2
            - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );
        
    scalar reScale = sqrt(EInc/magSqr(U));
    
    U *= reScale;
    
//     Info << "Energy after 1 = " << magSqr(U) << endl;
    
    U += uNew;
    
    measurePropertiesAfterControl(p, 0.0);
    
//     Info << "Energy after 2 = " << magSqr(U) << endl;
}

void dsmcAdiabaticWallRotationPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcAdiabaticWallRotationPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}



} // End namespace Foam

// ************************************************************************* //
