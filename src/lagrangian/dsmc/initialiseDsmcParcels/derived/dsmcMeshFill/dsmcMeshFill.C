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

#include "dsmcMeshFill.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcMeshFill, 0);

addToRunTimeSelectionTable(dsmcConfiguration, dsmcMeshFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMeshFill::dsmcMeshFill
(
    dsmcCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    dsmcConfiguration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMeshFill::~dsmcMeshFill()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcMeshFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const scalar translationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("translationalTemperature"))
    );
    
    const scalar rotationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("rotationalTemperature"))
    );
    
    const scalar vibrationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("vibrationalTemperature"))
    );
    
    const scalar electronicTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("electronicTemperature"))
    );

    const vector velocity(dsmcInitialiseDict_.lookup("velocity"));

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict_.subDict("numberDensities")
    );

    wordList molecules(numberDensitiesDict.toc());

    scalarList numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
    }

    forAll(mesh_.cells(), cellI)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            cellI
        );

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
                
                particlesRequired /= cloud_.nParticles(cellI);

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    nParticlesToInsert++;
                }

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = cloud_.equipartitionLinearVelocity
                    (
                        translationalTemperature,
                        cP.mass()
                    );

                    scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        rotationalTemperature,
                        cP.rotationalDegreesOfFreedom()
                    );
            
                    labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                    (
                        vibrationalTemperature,
                        cP.nVibrationalModes(),
                        typeId
                    );
                    
                    label ELevel = cloud_.equipartitionElectronicLevel
                    (
                        electronicTemperature,
                        cP.electronicDegeneracyList(),
                        cP.electronicEnergyList(),
                        typeId
                    );

                    U += velocity;
                    
                    label newParcel = -1;
                    
                    label classification = 0;
                    
                    const scalar& RWF = 
                        cloud_.coordSystem().recalculateRWF(cellI);

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

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const dsmcParcel::constantProperties& cP = cloud_.constProps
    (
        mostAbundantType
    );

    cloud_.sigmaTcRMax().primitiveFieldRef() = cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
    (
        translationalTemperature,
        cP.mass()
    );

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


} // End namespace Foam

// ************************************************************************* //
