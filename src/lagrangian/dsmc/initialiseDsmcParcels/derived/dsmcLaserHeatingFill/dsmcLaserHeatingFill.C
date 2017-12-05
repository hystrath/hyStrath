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

#include "dsmcLaserHeatingFill.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcLaserHeatingFill, 0);

addToRunTimeSelectionTable(dsmcConfiguration, dsmcLaserHeatingFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcLaserHeatingFill::dsmcLaserHeatingFill
(
    dsmcCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    dsmcConfiguration(cloud, dict),
    r0_(0.0),
    deltaT_(0.0),
    centrePoint_(vector::zero)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcLaserHeatingFill::~dsmcLaserHeatingFill()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcLaserHeatingFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

//     const scalar translationalTemperature
//     (
//         readScalar(dsmcInitialiseDict_.lookup("translationalTemperature"))
//     );
//     
//     const scalar rotationalTemperature
//     (
//         readScalar(dsmcInitialiseDict_.lookup("rotationalTemperature"))
//     );
//     
//     const scalar vibrationalTemperature
//     (
//         readScalar(dsmcInitialiseDict_.lookup("vibrationalTemperature"))
//     );
//     
//     const scalar electronicTemperature
//     (
//         readScalar(dsmcInitialiseDict_.lookup("electronicTemperature"))
//     );

    const vector velocity(dsmcInitialiseDict_.lookup("velocity"));
    
    r0_ = readScalar((dsmcInitialiseDict_.lookup("r0")));
    
    centrePoint_ = (dsmcInitialiseDict_.lookup("centrePoint"));
    
    deltaT_ = readScalar((dsmcInitialiseDict_.lookup("axialTemperature")));

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict_.subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
    }

    numberDensities /= cloud_.nParticle();

    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word regionName(dsmcInitialiseDict_.lookup("zoneName"));
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("dsmcLaserHeatingFill::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"dsmcInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    if (zone.size())
    {
        Info << "Lattice in zone: " << regionName << endl;

        forAll(zone, c)
        {
            const label& cellI = zone[c];
    
            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh_,
                cellI
            );
            
            const scalar& cellCentreI = mag(mesh_.cellCentres()[cellI] - centrePoint_);
            
            const scalar& cellCentreY = mesh_.cellCentres()[cellI].y();
            
            scalar cellTemperature = deltaT_*exp(-1.0*pow(cellCentreY/r0_,2.0));
            
            const scalar& cellCentreX = mesh_.cellCentres()[cellI].x();
            
            cellTemperature *= (0.0218 - cellCentreX)/0.0036;
    
            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];
    
                tetPointRef tet = cellTetIs.tet(mesh_);
    
                scalar tetVolume = tet.mag();
    
                forAll(molecules, i)
                {
                    const word& moleculeName(molecules[i]);
    
                    label typeId(findIndex(cloud_.typeIdList(), moleculeName));
    
                    if (typeId == -1)
                    {
                        FatalErrorIn("Foam::dsmcCloud<dsmcParcel>::initialise")
                            << "typeId " << moleculeName << "not defined." << nl
                            << abort(FatalError);
                    }
    
                    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);
    
                    scalar numberDensity = numberDensities[i];
    
                    // Calculate the number of particles required
                    scalar particlesRequired = numberDensity*tetVolume;
                    
                    const scalar& RWF = cloud_.getRWF_cell(cellI);
                    particlesRequired /= RWF;
    
                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);
    
                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                            > rndGen_.scalar01()
                    )
                    {
                        nParticlesToInsert++;
                    }
    
                    for (label pI = 0; pI < nParticlesToInsert; pI++)
                    {
                        point p = tet.randomPoint(rndGen_);
    
                        vector U = cloud_.equipartitionLinearVelocity
                        (
                            cellTemperature,
                            cP.mass()
                        );
    
                        scalar ERot = cloud_.equipartitionRotationalEnergy
                        (
                            cellTemperature,
                            cP.rotationalDegreesOfFreedom()
                        );
            
                        labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                        (
                            cellTemperature,
                            cP.vibrationalDegreesOfFreedom(),
                            typeId
                        );
                        
                        label ELevel = cloud_.equipartitionElectronicLevel
                        (
                            cellTemperature,
                            cP.degeneracyList(),
                            cP.electronicEnergyList(),
                            typeId
                        );

                        U += velocity;
                        
                        label newParcel = 0;
                        
                        label classification = 0;
                        
                        const scalar& RWF = cloud_.getRWF_cell(cellI);
                        
                        cloud_.addNewParcel
                        (
                            p,
                            U,
                            RWF,
                            ERot,
                            ELevel,
                            cellI,
                            cellTetIs.face(),
                            cellTetIs.tetPt(),
                            typeId,
                            newParcel,
                            classification,
                            vibLevel
                        );
                    }
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const dsmcParcel::constantProperties& cP = cloud_.constProps
    (
        mostAbundantType
    );


    forAll(zone, c)
    {
        const label& cellI = zone[c];
        
        const scalar& cellCentreI = mag(mesh_.cellCentres()[cellI] - centrePoint_);
            
        scalar cellTemperature = deltaT_*exp(-1.0*sqr(cellCentreI/r0_));

        cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] = cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
        (
            cellTemperature,
            cP.mass()
        );
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


} // End namespace Foam

// ************************************************************************* //
