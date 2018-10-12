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

#include "pdByField.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdByField, 0);

addToRunTimeSelectionTable(pdConfiguration, pdByField, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdByField::pdByField
(
    pdCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    pdConfiguration(cloud, dict),
    numberDensities_()
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdByField::~pdByField()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void pdByField::setInitialConfiguration()
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

    //- initialise numberDensities field. Outer list = species, inner = scalarField
    const List<word> molecules(pdInitialiseDict_.lookup("typeIds"));

    const List<word> densityFieldNames_(pdInitialiseDict_.lookup("densityFields"));

    const List<word> velocityFieldNames_(pdInitialiseDict_.lookup("velocityFields"));

    numberDensities_.setSize(molecules.size());
    particleVelocity_.setSize(molecules.size());

    forAll(numberDensities_, i)
    {
        numberDensities_[i].setSize(mesh_.nCells());
        particleVelocity_[i].setSize(mesh_.nCells());

        Info << "Building " << molecules[i] << " distribution from: " << endl;

        Info << "        density: " << densityFieldNames_[i] << endl;

        volScalarField rho
        (
            IOobject
            (
                densityFieldNames_[i],
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        if(i < velocityFieldNames_.size())
        {
            Info << "        velocity: " << velocityFieldNames_[i] << endl;
            volVectorField U
            (
                IOobject
                (
                    velocityFieldNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

            forAll(numberDensities_[i],cI)
            {
                numberDensities_[i][cI] = rho.internalField()[cI]/cloud_.nParticle();
                particleVelocity_[i][cI] = U.internalField()[cI];
            }
        }
        else
        {
            Info << "        velocity: no field provide, using uniform velocity " << velocity << endl;
            volVectorField U
            (
                IOobject
                (
                    "U",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector("0.0", dimLength/dimTime, velocity)
            );

            forAll(numberDensities_[i],cI)
            {
                numberDensities_[i][cI] = rho.internalField()[cI]/cloud_.nParticle();
                particleVelocity_[i][cI] = U.internalField()[cI];
            }
        }
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

                label typeId(findIndex(cloud_.typeIdList(), molecules[i]));

                if (typeId == -1)
                {
                    FatalErrorIn("Foam::pdCloud<pdParcel>::initialise")
                        << "typeId " << molecules[i] << "not defined." << nl
                        << abort(FatalError);
                }

                const pdParcel::constantProperties& cP = cloud_.constProps(typeId);

                scalar numberDensity = numberDensities_[i][cellI];

                // Calculate the number of particles required
                scalar particlesRequired = numberDensity*tetVolume;

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

                    scalar EVib = cloud_.equipartitionVibrationalEnergy
                    (
                        vibrationalTemperature,
                        cP.vibrationalDegreesOfFreedom(),
                        typeId
                    );

                    //initialising acceleration and force on particle as [0;0;0] at the moment
                    vector A = vector::zero;
                    scalar EPot = 0.0;

                    U += particleVelocity_[i][cellI];

                    label newParcel = 1;

                    label classification = 0;

                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        A,
                        EPot,
                        ERot,
                        EVib,
                        cellI,
                        cellTetIs.face(),
                        cellTetIs.tetPt(),
                        typeId,
                        newParcel,
                        classification
                    );
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)
    label mostAbundantType = -1;

    scalar temp = 0.0;
    forAll(numberDensities_,i)
    {
        if(findMax(numberDensities_[i]) > temp)
        {
            temp = findMax(numberDensities_[i]);
            mostAbundantType = i;
        }
    }


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
