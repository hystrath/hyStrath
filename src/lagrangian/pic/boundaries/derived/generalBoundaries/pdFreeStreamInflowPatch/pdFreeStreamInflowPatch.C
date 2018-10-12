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

#include "pdFreeStreamInflowPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(pdFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable(pdGeneralBoundary, pdFreeStreamInflowPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdFreeStreamInflowPatch::pdFreeStreamInflowPatch
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    numberDensities_(),
    accumulatedParcelsToInsert_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdFreeStreamInflowPatch::~pdFreeStreamInflowPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void pdFreeStreamInflowPatch::initialConfiguration()
{}

void pdFreeStreamInflowPatch::calculateProperties()
{

}

void pdFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();
    const scalar deltaT = mesh_.time().deltaTValue();

    scalar sqrtPi = sqrt(pi);


    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label& typeId = typeIds_[i];
        scalar mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    translationalTemperature_,
                    mass
                )
            );

             // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity_ & -sF/fA )/mostProbableSpeed;

            // From Bird eqn 4.22

            accumulatedParcelsToInsert_[i][f] +=
            (
                fA*numberDensities_[i]*deltaT*mostProbableSpeed
                *
                (
                   exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
            /(2.0*sqrtPi*cloud_.nParticle());
        }
    }


    labelField parcelsInserted(typeIds_.size(), 0);
    labelField parcelsToAdd(typeIds_.size(), 0);

    const vector& faceVelocity = velocity_;
    const scalar& faceTranslationalTemperature = translationalTemperature_;
    const scalar& faceRotationalTemperature = rotationalTemperature_;
    const scalar& faceVibrationalTemperature = vibrationalTemperature_;

    // insert pacels
    forAll(faces_, f)
    {
        const label& faceI = faces_[f];
        const label& cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

//         Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), 0.0);

        scalar previousCummulativeSum = 0.0;

        forAll(faceTets, triI)
        {
            const tetIndices& faceTetIs = faceTets[triI];

            cTriAFracs[triI] =
                faceTetIs.faceTri(mesh_).mag()/fA
                + previousCummulativeSum;

            previousCummulativeSum = cTriAFracs[triI];
        }

//         Force the last area fraction value to 1.0 to avoid any
//         rounding/non-flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

//         Normal unit vector *negative* so normal is pointing into the
//         domain
        vector n = sF;
        n /= -mag(n);

//         Wall tangential unit vector. Use the direction between the
//         face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

//         Other tangential unit vector.  Rescaling in case face is not
//         flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, m)
        {
            const label& typeId = typeIds_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];

            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                nI++;
            }

            faceAccumulator -= nI;
            parcelsToAdd[m] += nI;

            scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nI; i++)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5*
                    (
                        1.0
                        + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's PD3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if(abs(faceVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                                *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()*faceTranslationalTemperature/mass)
                    *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                    + (t1 & faceVelocity)*t1
                    + (t2 & faceVelocity)*t2
                    + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );

                scalar EVib = cloud_.equipartitionVibrationalEnergy
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).vibrationalDegreesOfFreedom(),
                    typeId
                );

                label newParcel = 1;

                /***********************************************************************/
                //initialising accerleration and force on particle as (0, 0, 0) at the moment
                vector A = vector::zero;
                scalar EPot = 0.0;
                /***********************************************************************/

                cloud_.addNewParcel
                (
                    p,
                    U,
                    /***********************************************************************/
                    A,
                    EPot,
                    /***********************************************************************/
                    ERot,
                    EVib,
                    cellI,
                    faces_[f],
                    faceTetIs.tetPt(),
                    typeId,
                    newParcel,
                    0
                );

                parcelsInserted[m] += 1.0;
            }
        }
    }


    if (Pstream::parRun())
    {
        forAll(parcelsInserted, m)
        {
            reduce(parcelsToAdd[m], sumOp<scalar>());
            reduce(parcelsInserted[m], sumOp<scalar>());

            Info<< "Specie: " << typeIds_[m]
                << ", target parcels to insert: " << parcelsToAdd[m]
                <<", inserted parcels: " << parcelsInserted[m]
                << endl;
        }
    }
}

void pdFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{

}

void pdFreeStreamInflowPatch::controlParcelsAfterCollisions()
{
}

void pdFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void pdFreeStreamInflowPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

//     setProperties();
}



void pdFreeStreamInflowPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    translationalTemperature_ = readScalar(propsDict_.lookup("translationalTemperature"));
    rotationalTemperature_ = readScalar(propsDict_.lookup("rotationalTemperature"));
    vibrationalTemperature_ = readScalar(propsDict_.lookup("vibrationalTemperature"));

    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("pdFreeStreamInflowPatch::pdFreeStreamInflowPatch()")
            << "Cannot have zero typeIds being inserd." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("pdFreeStreamInflowPatch::pdFreeStreamInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // read in the mass density per specie

    const dictionary& numberDensitiesDict
    (
        propsDict_.subDict("numberDensities")
    );

    numberDensities_.clear();

    numberDensities_.setSize(typeIds_.size(), 0.0);

    forAll(numberDensities_, i)
    {
        numberDensities_[i] = readScalar
        (
            numberDensitiesDict.lookup(moleculesReduced[i])
        );
    }

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}



} // End namespace Foam

// ************************************************************************* //
