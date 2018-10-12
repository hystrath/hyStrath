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

#include "pdZoneFill.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdZoneFill, 0);

addToRunTimeSelectionTable(pdConfiguration, pdZoneFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdZoneFill::pdZoneFill
(
    pdCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    pdConfiguration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdZoneFill::~pdZoneFill()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void pdZoneFill::setInitialConfiguration()
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

    const dictionary& numberDensitiesDict
    (
        pdInitialiseDict_.subDict("numberDensities")
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
    const word regionName(pdInitialiseDict_.lookup("zoneName"));
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("pdZoneFill::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"pdInitialiseDict"
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
                        FatalErrorIn("Foam::pdCloud<pdParcel>::initialise")
                            << "typeId " << moleculeName << "not defined." << nl
                            << abort(FatalError);
                    }

                    const pdParcel::constantProperties& cP = cloud_.constProps(typeId);

                    scalar numberDensity = numberDensities[i];

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

                        //initialising acceleration and force on particle as [0;0;0] at the moment
                        vector A = vector::zero;
                        scalar EPot = 0.0;

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
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const pdParcel::constantProperties& cP = cloud_.constProps
    (
        mostAbundantType
    );


    forAll(zone, c)
    {
        const label& cellI = zone[c];

        cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] = cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
        (
            translationalTemperature,
            cP.mass()
        );
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();

}


} // End namespace Foam

// ************************************************************************* //
