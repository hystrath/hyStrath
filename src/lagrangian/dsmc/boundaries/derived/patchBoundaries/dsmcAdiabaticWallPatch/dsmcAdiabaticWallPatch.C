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

#include "dsmcAdiabaticWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcAdiabaticWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary, 
    dsmcAdiabaticWallPatch, 
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcAdiabaticWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("initialTemperature"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcAdiabaticWallPatch::dsmcAdiabaticWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
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

dsmcAdiabaticWallPatch::~dsmcAdiabaticWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAdiabaticWallPatch::initialConfiguration()
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


void dsmcAdiabaticWallPatch::calculateProperties()
{            
    const scalar deltaT = mesh_.time().deltaTValue(); //TODO cloud_.deltaTValue(p.cell());
    
    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();
    
    forAll(cells_, c)
    {
        const label celli = cells_[c];
        
        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[celli];
                  
        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            
            const dsmcParcel::constantProperties& constProp 
                = cloud_.constProps(p->typeId());
                                
            const scalar& mass = constProp.mass()*cloud_.nParticles(celli);
                        
            totalMass_ += mass;
            totalMols_ += 1.0;
            mcc_ += mass*mag(p->U())*mag(p->U());
            incidentMomentum_ += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*p->U();
            incidentMass_ += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass();
                                
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
    
        const scalar gasTemperature = (1.0/(3.0*physicoChemical::k.value()))
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


void dsmcAdiabaticWallPatch::controlParticle
(
    dsmcParcel& p, 
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);
    
    const scalar mass = cloud_.constProps(p.typeId()).mass();
    
    vector& U = p.U();
    
    const scalar& magUInitial = mag(U);

    //- Wall unit normal vector and wall unit tangential vectors
    vector nw, tw1, tw2 = vector::zero;

    dsmcPatchBoundary::calculateWallUnitVectors(p, nw, tw1, tw2);

    Random& rndGen(cloud_.rndGen()); 

    const scalar Rf1 = rndGen.sample01<scalar>();
    const scalar Rf2 = rndGen.sample01<scalar>();
    
    U =  -magUInitial*nw*sqrt(Rf1)
        + magUInitial*sqrt(1.0 - Rf1)*tw1*cos(twoPi*Rf2)
        + magUInitial*sqrt(1.0 - Rf1)*tw2*sin(twoPi*Rf2);
        
    U += velocity_;

    measurePropertiesAfterControl(p, 0.0);
}

void dsmcAdiabaticWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcAdiabaticWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
