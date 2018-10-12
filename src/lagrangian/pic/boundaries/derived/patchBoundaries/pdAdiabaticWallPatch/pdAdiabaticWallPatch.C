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

#include "pdAdiabaticWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdAdiabaticWallPatch, 0);

addToRunTimeSelectionTable(pdPatchBoundary, pdAdiabaticWallPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdAdiabaticWallPatch::pdAdiabaticWallPatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    incidentMomentum_(vector::zero),
    gasVelocity_(vector::zero),
    totalVolume_(0.0),
    totalMass_(0.0),
    incidentMass_(0.0),
    wallTemperature_(0.0),
    mols_(0.0),
    totalMols_(0.0),
    cv_(0.0),
    mcc_(0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
    
    wallTemperature_ = temperature_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdAdiabaticWallPatch::~pdAdiabaticWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void pdAdiabaticWallPatch::initialConfiguration()
{
    forAll(cells_, c)
    {
        const label& cellI = cells_[c];
        totalVolume_ += mesh_.cellVolumes()[cellI];
    }
    
    if (Pstream::parRun())
    {
        reduce(totalVolume_, sumOp<scalar>());
    }
}

void pdAdiabaticWallPatch::calculateProperties()
{            
    const scalar deltaT = mesh_.time().deltaTValue();
    
    const List<DynamicList<pdParcel*> >& cellOccupancy = cloud_.cellOccupancy();
    
    forAll(cells_, c)
    {
        const List<pdParcel*>& parcelsInCell = cellOccupancy[cells_[c]];
                  
        forAll(parcelsInCell, pIC)
        {
            pdParcel* p = parcelsInCell[pIC];
            
            const pdParcel::constantProperties& constProp 
                                = cloud_.constProps(p->typeId());
                                
            const scalar& mass = constProp.mass()*cloud_.nParticle();
                        
            totalMass_ += mass;
            totalMols_ += 1.0;
            mcc_ += mass*mag(p->U())*mag(p->U());
            incidentMomentum_ += cloud_.nParticle()*cloud_.constProps(p->typeId()).mass()*p->U();
            incidentMass_ += cloud_.nParticle()*cloud_.constProps(p->typeId()).mass();
                                
            //check if particle is heading towards the surface
            vector nw = p->normal();
            nw /= mag(nw);
                                
            // Normal velocity magnitude
            scalar U_dot_nw = p->U() & nw;
            
//             // Wall tangential velocity (flow direction)
//             vector Ut = p->U() - U_dot_nw*nw;
//             
//             // Wall tangential unit vector
//             vector tw1 = Ut/mag(Ut);
// 
//             // Other tangential unit vector
//             vector tw2 = nw^tw1;
            
            if(U_dot_nw < VSMALL)
            {                                      
                mols_ += 1.0;
                cv_ += magSqr(p->U() - gasVelocity_)*(U_dot_nw);
                
            }
        }
    }
    
    vector incidentMomentum = incidentMomentum_;
    scalar totalMass = totalMass_;
    scalar incidentMass = incidentMass_;
    scalar mols = mols_;
    scalar totalMols = totalMols_;
    scalar cv = cv_;
    scalar mcc = mcc_;
    
    if(Pstream::parRun())
    {
        reduce(incidentMomentum, sumOp<vector>());
        reduce(totalMass, sumOp<scalar>());
        reduce(incidentMass, sumOp<scalar>());
        reduce(mols, sumOp<scalar>());
        reduce(totalMols, sumOp<scalar>());
        reduce(cv, sumOp<scalar>());
        reduce(mcc, sumOp<scalar>());
    }
    
    if(mols > VSMALL)
    {
        gasVelocity_ = incidentMomentum/incidentMass;
    
        scalar gasTemperature = (1.0/(3.0*physicoChemical::k.value()))
                                *(
                                    ((mcc/(totalMols*cloud_.nParticle())))
                                    - (
                                        (totalMass/(totalMols*cloud_.nParticle())
                                        )*mag(gasVelocity_)*mag(gasVelocity_))
                                );
        
        Info << "gasVelocity = " << gasVelocity_ << endl;
        Info << "gasTemperature = " << gasTemperature << endl;
        
        scalar phi = -pow((2.0*1.3806e-23*gasTemperature/66.3e-27), 1.5)/sqrt(pi);
        scalar phi2 = 0.5*(cv/mols);
        
        Info << "phi2 = " << phi2 << endl;
        Info << "Equilibrium phi = " << phi << endl;
        
        Info << "Ratio = " << phi2 / (-pow((2.0*1.3806e-23*gasTemperature/66.3e-27), 1.5)/sqrt(pi)) << endl;
        
        vector Uslip = velocity_ - gasVelocity_;
        
        scalar err = GREAT;
        
        //newtonRaphson method
        do
        {
            scalar averageMass = totalMass/(totalMols*cloud_.nParticle());
            
            Info << "averageMass = " << averageMass << endl;
            
            scalar f = 1.0 + 
                    (averageMass*magSqr(Uslip))/(4.0*physicoChemical::k.value()*wallTemperature_) +
                    phi2*sqrt(pi)*pow(averageMass/(2.0*physicoChemical::k.value()*wallTemperature_),1.5);
                    
//             scalar df = (-averageMass/(4.0*physicoChemical::k.value()*sqr(wallTemperature_)))
//                     *(magSqr(Uslip) + 3.0*phi2*sqrt((pi*averageMass)/(2.0*physicoChemical::k.value()*wallTemperature_)));
                    
            scalar df = - (averageMass/(8.0*physicoChemical::k.value()*sqr(wallTemperature_)))*
                    (3.0*sqrt(twoPi)*phi2*sqrt(averageMass/(physicoChemical::k.value()*wallTemperature_)) + 2.0*magSqr(Uslip));
            
            Info << "f = " << f << endl;
            Info << "df = " << df << endl;
                
            scalar newTemperature = wallTemperature_ - f/df;
            
            err = fabs((wallTemperature_ - newTemperature)/*/wallTemperature_*/);
            
            Info << "err = " << err << endl;
            Info << "newTemperature = " << newTemperature << endl;
            
            wallTemperature_ = newTemperature;
                
        } while(wallTemperature_ > VSMALL && err > 0.1);
        
        Info << "newWallTemperature = " << wallTemperature_ << endl;
    }
    
    
    
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {        
//         incidentMomentum_ = vector::zero;
//         totalMass_ = 0.0;
//         incidentMass_ = 0.0;
//         mols_ = 0.0;
//         totalMols_ = 0.0;
//         cv_ = 0.0;
//         mcc_ = 0.0;  
    }
}

void pdAdiabaticWallPatch::controlParticle(pdParcel& p, pdParcel::trackingData& td)
{
    scalar mass = cloud_.constProps(p.typeId()).mass();
    
    vector& U = p.U();
    
    scalar magUInitial = mag(U);
    
    measurePropertiesBeforeControl(p);

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen()); 

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
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;
    
    scalar Rf1 = rndGen.sample01<scalar>();
    scalar Rf2 = rndGen.sample01<scalar>();
    
    U = -magUInitial*nw*sqrt(Rf1)
        + magUInitial*sqrt(1.0 - Rf1)*tw1*cos(twoPi*Rf2)
        + magUInitial*sqrt(1.0 - Rf1)*tw2*sin(twoPi*Rf2);

//     U =
//         sqrt(physicoChemical::k.value()*wallTemperature_/mass)
//         *(
//             rndGen.GaussNormal()*tw1
//             + rndGen.GaussNormal()*tw2
//             - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
//         );
        
    measurePropertiesAfterControl(p);
    
    U += velocity_;
    
//     incidentMomentum_ += cloud_.nParticle()*cloud_.constProps(p.typeId()).mass()*p.U();
//     incidentMass_ += cloud_.nParticle()*cloud_.constProps(p.typeId()).mass();
}

void pdAdiabaticWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{	
}


void pdAdiabaticWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();

}

void pdAdiabaticWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("initialTemperature"));
}

} // End namespace Foam

// ************************************************************************* //
