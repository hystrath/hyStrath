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

#include "pdInitialiseOnLine.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdInitialiseOnLine, 0);

addToRunTimeSelectionTable(pdConfiguration, pdInitialiseOnLine, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdInitialiseOnLine::pdInitialiseOnLine
(
    pdCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    pdConfiguration(cloud, dict),
    startPosition_(vector::zero),
    endPosition_(vector::zero)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdInitialiseOnLine::~pdInitialiseOnLine()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void pdInitialiseOnLine::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const scalar translationalTemperature
    (
        readScalar(pdInitialiseDict_.lookup("translationalTemperature"))
    );

    const scalar rotationalTemperature
    (
        readScalar(pdInitialiseDict_.lookup("rotationalTemperature"))
    );

    const scalar vibrationalTemperature
    (
        readScalar(pdInitialiseDict_.lookup("vibrationalTemperature"))
    );

    const vector velocity(pdInitialiseDict_.lookup("velocity"));

    startPosition_ = (pdInitialiseDict_.lookup("startPosition"));

    endPosition_ = (pdInitialiseDict_.lookup("endPosition"));

    scalar lineLength = mag(endPosition_ - startPosition_);

    vector unitVector = (endPosition_ - startPosition_)/mag(endPosition_ - startPosition_);

    const dictionary& numberDensitiesDict
    (
        pdInitialiseDict_.subDict("numberOfPDParticlesOnLine")
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

//     scalar nParticlesToInsert = 0.0;
//
//     forAll(molecules, i)
//     {
//         nParticlesToInsert += numberDensities[i];
//     }
//
//     scalar deltaPosition = mag(startPosition_ - endPosition_)/nParticlesToInsert;
//
//     vectorField positions(nParticlesToInsert, vector::zero);
//
//     forAll(positions, i)
//     {
//         positions[i] = startPosition_ + (0.5 + scalar(i))*deltaPosition*unitVector;
//     }

//     numberDensities /= cloud_.nParticle();

    forAll(molecules, i)
    {
        label nParticlesToInsert = label(numberDensities[i]);

        const word& moleculeName(molecules[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if (typeId == -1)
        {
            FatalErrorIn("Foam::pdCloud<pdParcel>::initialise")
                << "typeId " << moleculeName << "not defined." << nl
                << abort(FatalError);
        }

            const pdParcel::constantProperties& cP = cloud_.constProps(typeId);

        for (label pI = 0; pI < nParticlesToInsert; pI++)
        {
            point p = startPosition_ + (rndGen_.sample01<scalar>()*lineLength*unitVector);

            vector& findCellIndex = p;

            label cellI = mesh_.findCell(findCellIndex);

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

            scalar EVib = cloud_.equipartitionVibrationalEnergy
            (
                vibrationalTemperature,
                cP.vibrationalDegreesOfFreedom(),
                typeId
            );

            U += velocity;

            label newParcel = 1;

            label classification = 0;

            label tetFace = -1;
            label tetPt = -1;

            mesh_.findTetFacePt
            (
                cellI,
                p,
                tetFace,
                tetPt
            );

            //initialising acceleration and force on particle as [0;0;0] at the moment
            vector A = vector::zero;
            scalar EPot = 0.0;

            cloud_.addNewParcel
            (
                p,
                U,
                A,
                EPot,
                ERot,
                EVib,
                cellI,
                tetFace,
                tetPt,
                typeId,
                newParcel,
                classification
            );
        }
    }
//         List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
//         (
//             mesh_,
//             cellI
//         );
//
//         forAll(cellTets, tetI)
//         {
//             const tetIndices& cellTetIs = cellTets[tetI];
//
//             tetPointRef tet = cellTetIs.tet(mesh_);
//
//             scalar tetVolume = tet.mag();
//
//             forAll(molecules, i)
//             {
//                 const word& moleculeName(molecules[i]);
//
//                 label typeId(findIndex(cloud_.typeIdList(), moleculeName));
//
//                 if (typeId == -1)
//                 {
//                     FatalErrorIn("Foam::pdCloud<pdParcel>::initialise")
//                         << "typeId " << moleculeName << "not defined." << nl
//                         << abort(FatalError);
//                 }
//
//                 const pdParcel::constantProperties& cP = cloud_.constProps(typeId);
//
//                 scalar numberDensity = numberDensities[i];
//
//                 // Calculate the number of particles required
//                 scalar particlesRequired = numberDensity*tetVolume;
//
//                 // Only integer numbers of particles can be inserted
//                 label nParticlesToInsert = label(particlesRequired);
//
//                 // Add another particle with a probability proportional to the
//                 // remainder of taking the integer part of particlesRequired
//                 if
//                 (
//                     (particlesRequired - nParticlesToInsert)
//                   > rndGen_.scalar01()
//                 )
//                 {
//                     nParticlesToInsert++;
//                 }
//
//                 for (label pI = 0; pI < nParticlesToInsert; pI++)
//                 {
//                     point p = tet.randomPoint(rndGen_);
//
//                     vector U = cloud_.equipartitionLinearVelocity
//                     (
//                         translationalTemperature,
//                         cP.mass()
//                     );
//
//                     scalar ERot = cloud_.equipartitionRotationalEnergy
//                     (
//                         rotationalTemperature,
//                         cP.rotationalDegreesOfFreedom()
//                     );
//
//                     scalar EVib = cloud_.equipartitionVibrationalEnergy
//                     (
//                         vibrationalTemperature,
//                         cP.vibrationalDegreesOfFreedom(),
//                         typeId
//                     );
//
//                     U += velocity;
//
//                     label newParcel = 0;
//
//                     label classification = 0;
//
//                     cloud_.addNewParcel
//                     (
//                         p,
//                         U,
//                         ERot,
//                         EVib,
//                         cellI,
//                         cellTetIs.face(),
//                         cellTetIs.tetPt(),
//                         typeId,
//                         newParcel,
//                         classification
//                     );
//                 }
//             }
//         }
//     }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const pdParcel::constantProperties& cP = cloud_.constProps
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
