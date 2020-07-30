/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2020 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "dsmcMassFlowRateInlet.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcMassFlowRateInlet, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcMassFlowRateInlet, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcMassFlowRateInlet::dsmcMassFlowRateInlet
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    moleFractions_(),
    inletVelocity_(faces_.size(), vector::zero),
    previousInletVelocity_(faces_.size(), vector::zero),
    momentum_(faces_.size(), vector::zero),
    n_(),
    accumulatedParcelsToInsert_(),
    mass_(faces_.size(), 0.0),
    massFlowRate_(),
    inletTemperature_(),
//     theta_(),
    initialVelocity_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    //sorts issues with the velocity pointing out of the mesh
    inletVelocity_ = initialVelocity_;
    previousInletVelocity_ = initialVelocity_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcMassFlowRateInlet::~dsmcMassFlowRateInlet()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcMassFlowRateInlet::initialConfiguration()
{}

void dsmcMassFlowRateInlet::calculateProperties()
{

}

void dsmcMassFlowRateInlet::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    //loop over all species
    forAll(accumulatedParcelsToInsert_, iD)   // I Added.
    {
        // loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            vector faceVelocity = inletVelocity_[f];
            scalar faceTemperature = inletTemperature_;
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

            //Cumulative triangle area fractions
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

            //Force the last area fraction value to 1.0 to avoid any
            // rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //  Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            /* -----------------------------------------------------------------------------*/

            //generate Poisson distributed random number of particles to insert
            //see Tysanner & Garcia, Int. J. Numer. Meth. Fluids 2050; 00:1–12
            //this eliminates non-equilibrium behaviour that does not exist in the corresponding physical system

            label k = 0;

            const scalar target= exp(-accumulatedParcelsToInsert_[iD][f]);

            scalar p = rndGen.sample01<scalar>();

            while (p > target)
            {
                p *= rndGen.sample01<scalar>();
                k += 1;
            }

            label nParcelsToInsert = k;

            /* -----------------------------------------------------------------------------*/

//             label nParcelsToInsert = label(accumulatedParcelsToInsert_[iD][f]);
//
//             if ((nParcelsToInsert - accumulatedParcelsToInsert_[iD][f]) > rndGen.sample01<scalar>())
//             {
//                 nParcelsToInsert++;
//             }

            accumulatedParcelsToInsert_[iD][f] -= nParcelsToInsert; //remainder has been set

            const label& typeId = typeIds_[iD];
            scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; i++)
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
                        faceTemperature,
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

                // Equivalent to the QA value in Bird's DSMC3.FOR
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
                    sqrt(physicoChemical::k.value()*faceTemperature/mass)
                    *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                    + (t1 & faceVelocity)*t1
                    + (t2 & faceVelocity)*t2
                    + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).nVibrationalModes(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).electronicDegeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList()
                );

                label newParcel = patchId();

                const scalar& RWF = cloud_.coordSystem().RWF(cellI);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    faces_[f],
                    faceTetIs.tetPt(),
                    typeId,
                    newParcel,
                    0,
                    vibLevel
                );
            }
        }
    }
    previousInletVelocity_ = inletVelocity_;
}

void dsmcMassFlowRateInlet::controlParcelsBeforeCollisions()
{
    // insert pacels after move, but before collisions - REMINDER: NEW PARCELS MUST BE ADDED TO CELL OCCUPANCY MANUALLY


}

void dsmcMassFlowRateInlet::controlParcelsAfterCollisions()
{
    const scalar sqrtPi = sqrt(pi);

//     vectorField momentum(faces_.size(), vector::zero);
    vectorField newInletVelocity(faces_.size(), vector::zero);
//     scalarField mass(faces_.size(), scalar(0.0));

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];

        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[celli];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];

            momentum_[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*p->U();
            mass_[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass();
        }

//         newInletVelocity[c] = momentum[c]/mass[c];

//         inletVelocity_[c] = theta_*newInletVelocity[c] + (1.0 - theta_)*previousInletVelocity_[c];

        inletVelocity_[c] = momentum_[c]/mass_[c];

        const vector& sF = mesh_.faceAreas()[faces_[c]];
        const scalar fA = mag(sF);

        if((inletVelocity_[c] & -sF/fA) < 0)
        {
            inletVelocity_[c] = previousInletVelocity_[c];
        }

//         if( (inletVelocity_[c] & -sF/fA) > (100.0*(initialVelocity_ & -sF/fA)) )
//         {
//             inletVelocity_[c] = previousInletVelocity_[c];
//         }
    }

    if(faces_.size() > VSMALL)
    {
        Pout << "dsmcMassFlowRateInlet inlet velocity = " << inletVelocity_[(faces_.size()/2)] << endl;
    }

    scalarField massFractions(typeIds_.size(), 0.0);
    scalar totalMass = 0.0;

    forAll(massFractions, iD)
    {
        const label& typeId = typeIds_[iD];

        scalar mass = cloud_.constProps(typeId).mass();

        totalMass += mass*moleFractions_[iD];
    }

    forAll(massFractions, iD)
    {
        const label& typeId = typeIds_[iD];

        scalar mass = cloud_.constProps(typeId).mass();

        massFractions[iD] = moleFractions_[iD]*(mass/totalMass);
    }

    Info << "dsmcMassFlowRateInlet massFractions = " << massFractions << endl;

    // compute number of parcels to insert
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label& typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT = cloud_.deltaTValue(mesh_.boundaryMesh()[patchId_].faceCells()[faceI]);

            scalar mass = cloud_.constProps(typeId).mass();

            n_[iD][f] = (massFractions[iD]*massFlowRate_)/
                    ((inletVelocity_[f] & -sF/fA)*patchSurfaceArea_*mass);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    inletTemperature_,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (inletVelocity_[f] & -sF/fA )/mostProbableSpeed;

            // From Bird eqn 4.22
            accumulatedParcelsToInsert_[iD][f] +=
//                 moleFractions_[iD]*
                (
                    fA*n_[iD][f]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
        }
    }
}

void dsmcMassFlowRateInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcMassFlowRateInlet::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcMassFlowRateInlet::setProperties()
{
    massFlowRate_ = readScalar(propsDict_.lookup("massFlowRate"));

    inletTemperature_ = readScalar(propsDict_.lookup("inletTemperature"));

    initialVelocity_ = propsDict_.lookup("initialVelocity");

//     theta_ = readScalar(propsDict_.lookup("theta"));

//     if(0.0 > theta_ || theta_ > 1.0)
//     {
//         FatalErrorIn("dsmcLiouFangPressureInlet::dsmcLiouFangPressureInlet()")
//             << "Theta must be a value between 0 and 1 " << nl << "in: "
//             << mesh_.time().system()/"boundariesDict"
//             << exit(FatalError);
//     }

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcMassFlowRateInlet::dsmcMassFlowRateInlet()")
            << "Cannot have zero typeIds being inserted." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, iD)
    {
        const word& moleculeName(molecules[iD]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, iD)
    {
        const word& moleculeName(moleculesReduced[iD]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcMassFlowRateInlet::setProperties()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[iD] = typeId;
    }

    // read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, iD)
    {
        moleFractions_[iD] = readScalar
        (
            moleFractionsDict.lookup(moleculesReduced[iD])
        );
    }

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }


    n_.setSize(typeIds_.size());

    forAll(n_, m)
    {
        n_[m].setSize(nFaces_, 0.0);
    }

}


void dsmcMassFlowRateInlet::setNewBoundaryFields()
{

}



} // End namespace Foam

// ************************************************************************* //
