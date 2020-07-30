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


\*--------------------------------------------------------------------------*/

#include "dsmcLiouFangPressureOutletSpecifiedMolarFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcLiouFangPressureOutletSpecifiedMolarFraction, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcLiouFangPressureOutletSpecifiedMolarFraction, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcLiouFangPressureOutletSpecifiedMolarFraction::dsmcLiouFangPressureOutletSpecifiedMolarFraction
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    outletPressure_(),
    nTimeSteps_(scalar(0.0)),
    typeIds_(),
    UMean_(faces_.size(), vector::zero),
    outletVelocity_(faces_.size(), vector::zero),
    UCollected_(faces_.size(), vector::zero),
    accumulatedParcelsToInsert_(),
    totalVibrationalEnergy_(),
    vibT_(),
    vDof_(),
    nTotalParcelsSpecies_(),
    moleFractions_(),
    cellVolume_(faces_.size(), scalar(0.0)),
    outletNumberDensity_(faces_.size(), scalar(0.0)),
    outletMassDensity_(faces_.size(), scalar(0.0)),
    outletTemperature_(faces_.size(), scalar(300.0)),
    totalMomentum_(faces_.size(), vector::zero),
    totalMass_(faces_.size(), scalar(0.0)),
    totalMassDensity_(faces_.size(), scalar(0.0)),
    totalNumberDensity_(faces_.size(), scalar(0.0)),
    totalRotationalEnergy_(faces_.size(), scalar(0.0)),
    totalRotationalDof_(faces_.size(), scalar(0.0)),
    nTotalParcels_(faces_.size(), scalar(0.0)),
    nTotalParcelsInt_(faces_.size(), scalar(0.0)),
    mcc_(faces_.size(), scalar(0.0))
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];  //get volume of each boundary cell
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcLiouFangPressureOutletSpecifiedMolarFraction::~dsmcLiouFangPressureOutletSpecifiedMolarFraction()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcLiouFangPressureOutletSpecifiedMolarFraction::initialConfiguration()
{}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::calculateProperties()
{

}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    //loop over all species
    forAll(accumulatedParcelsToInsert_, iD)   // I Added.
    {
        // loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            vector faceVelocity = outletVelocity_[f];
            scalar faceTemperature = outletTemperature_[f];
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
            //rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //Wall tangential unit vector. Use the direction between the
            //face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            label nParcelsToInsert = label(accumulatedParcelsToInsert_[iD][f]);

            if ((nParcelsToInsert - accumulatedParcelsToInsert_[iD][f]) > rndGen.sample01<scalar>())
            {
                nParcelsToInsert++;
            }

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

                label newParcel = patchId();

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).electronicDegeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList()
                );

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
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsBeforeCollisions()
{

}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    const scalar sqrtPi = sqrt(pi);

    scalar molecularMass = 0.0;
    scalar molarcontantPressureSpecificHeat = 0.0;
    scalar molarcontantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label& typeId = typeIds_[iD];

        molecularMass += cloud_.constProps(typeId).mass()*moleFractions_[iD];
        molarcontantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[iD];
        molarcontantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[iD];
    }

    // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma = molarcontantPressureSpecificHeat/molarcontantVolumeSpecificHeat;

    // calculate properties in cells attached to each boundary face

    vectorField momentum(faces_.size(), vector::zero);
    vectorField UCollected(faces_.size(), vector::zero);
    scalarField mass(faces_.size(), scalar(0.0));
    scalarField nParcels(faces_.size(), scalar(0.0));
    scalarField nParcelsInt(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanX(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanY(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanZ(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrX(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrY(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrZ(faces_.size(), scalar(0.0));
    scalarField rotationalEnergy(faces_.size(), scalar(0.0));
    scalarField rotationalDof(faces_.size(), scalar(0.0));
    scalarField mcc(faces_.size(), scalar(0.0));
    scalarField vDoF(faces_.size(), scalar(0.0));
    scalarField translationalTemperature(faces_.size(), scalar(0.0));
    scalarField rotationalTemperature(faces_.size(), scalar(0.0));
    scalarField vibrationalTemperature(faces_.size(), scalar(0.0));
    scalarField overallTemperature(faces_.size(), scalar(0.0));
    scalarField numberDensity(faces_.size(), scalar(0.0));
    scalarField massDensity(faces_.size(), scalar(0.0));
    scalarField pressure(faces_.size(), scalar(0.0));
    scalarField speedOfSound(faces_.size(), scalar(0.0));
    scalarField velocityCorrection(faces_.size(), scalar(0.0));
    scalarField massDensityCorrection(faces_.size(), scalar(0.0));

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];

        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];

            momentum[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*p->U();
            mass[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass();
            mcc[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*mag(p->U())*mag(p->U());
            nParcels[c] += 1.0;
            UCollected[c] += p->U();
            rotationalEnergy[c] += p->ERot();
            rotationalDof[c] += cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom();
//             velSqrMeanX[c] += sqr(p->U().x());
//             velSqrMeanY[c] += sqr(p->U().y());
//             velSqrMeanZ[c] += sqr(p->U().z());
//             velMeanSqrX[c] += p->U().x();
//             velMeanSqrY[c] += p->U().y();
//             velMeanSqrZ[c] += p->U().z();

            if(cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom() > VSMALL)
            {
                nParcelsInt[c] += 1.0;
            }
        }

        nTotalParcels_[c] += nParcels[c];

        nTotalParcelsInt_[c] += nParcelsInt[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        mcc_[c] += mcc[c];

        UCollected_[c] += UCollected[c];

        totalRotationalEnergy_[c] += rotationalEnergy[c];

        totalRotationalDof_[c] += rotationalDof[c];

        massDensity[c] = totalMass_[c] / (cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if(nTotalParcels_[c] > 1)
        {

//             velSqrMeanX_[c] += velSqrMeanX[c];
//             velSqrMeanY_[c] += velSqrMeanY[c];
//             velSqrMeanZ_[c] += velSqrMeanZ[c];
//             velMeanSqrX_[c] += velMeanSqrX[c];
//             velMeanSqrY_[c] += velMeanSqrY[c];
//             velMeanSqrZ_[c] += velMeanSqrZ[c];
//
//             totalVelSqrMeanX_[c] = velSqrMeanX_[c]/nTotalParcels_[c];
//             totalVelSqrMeanY_[c] = velSqrMeanY_[c]/nTotalParcels_[c];
//             totalVelSqrMeanZ_[c] = velSqrMeanZ_[c]/nTotalParcels_[c];
//             totalVelMeanSqrX_[c] = sqr(velMeanSqrX_[c]/nTotalParcels_[c]);
//             totalVelMeanSqrY_[c] = sqr(velMeanSqrY_[c]/nTotalParcels_[c]);
//             totalVelMeanSqrZ_[c] = sqr(velMeanSqrZ_[c]/nTotalParcels_[c]);
//
//             translationalTemperature[c] = (0.5*molecularMass)*(2.0/(3.0*physicoChemical::k.value()))
//                         *(
//                                 totalVelSqrMeanX_[c]+totalVelSqrMeanY_[c]+totalVelSqrMeanZ_[c]
//                                 -totalVelMeanSqrX_[c]-totalVelMeanSqrY_[c]-totalVelMeanSqrZ_[c]
//                             );

            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            translationalTemperature[c] = (1.0/(3.0*physicoChemical::k.value()))
                  *(
                      ((mcc_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))))
                      - (
                          (mass[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))
                        )*mag(UMean_[c])*mag(UMean_[c]))
                  );

//             Info << "translationalTemperature[c] = "  << translationalTemperature << endl;


            if(translationalTemperature[c] < VSMALL)
            {
                translationalTemperature[c] = 300.00;
            }

//             if(totalRotationalDof_[c] > VSMALL)
//             {
//                 rotationalTemperature[c] = (2.0/physicoChemical::k.value())*(totalRotationalEnergy_[c]/totalRotationalDof_[c]);
//             }
//             else
//             {
//                 rotationalTemperature[c] = 0.0;
//             }

//             forAll(totalVibrationalEnergy_, iD)
//             {
//                 const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];
//
//                 forAll(parcelsInCell, pIC)
//                 {
//                     dsmcParcel* p = parcelsInCell[pIC];
//
//                     if(p->typeId() == typeIds_[iD])
//                     {
//                         totalVibrationalEnergy_[iD][c] += p->vibLevel()*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV();
//                         nTotalParcelsSpecies_[iD][c] += 1.0;
//                     }
//                 }
//             }

//             forAll(totalVibrationalEnergy_, iD)
//             {
//                 if(totalVibrationalEnergy_[iD][c] > VSMALL)
//                 {
//                     const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
//
//                     scalar vibrationalEMean = (totalVibrationalEnergy_[iD][c]/nTotalParcelsSpecies_[iD][c]);
//                     scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//
//                     vibT_[iD][c] = thetaV / log(1.0 + (1.0/iMean));
//                     vDof_[iD][c] = (2.0*thetaV/vibT_[iD][c]) / (exp(thetaV/vibT_[iD][c]) - 1.0);
//
//                     scalar fraction = nTotalParcelsSpecies_[iD][c]/nTotalParcelsInt_[c];
//
//                     vDoF[c] += fraction*vDof_[iD][c];
//
//                     vibrationalTemperature[c] += fraction*vibT_[iD][c];
//                 }
//             }

            scalar nRotDof = 0.0;

            if(nTotalParcels_[c] > VSMALL)
            {
                nRotDof = totalRotationalDof_[c]/nTotalParcels_[c];
            }


//             overallTemperature[c] = (
//                                         (3.0*translationalTemperature[c])
//                                         + (nRotDof*rotationalTemperature[c])
//                                         + (vDoF[c]*vibrationalTemperature[c])
//                                     )
//                         / (3.0 + nRotDof + vDoF[c]);

            pressure[c] = numberDensity[c]*physicoChemical::k.value()*translationalTemperature[c];

            speedOfSound[c] = sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n /= mag(n);

            scalar faceNormalVelocity = (n & outletVelocity_[c]);

            if(nTimeSteps_ > 100)
            {
//             if(faceNormalVelocity < VSMALL)
//             {
                massDensityCorrection[c] = (outletPressure_ - pressure[c]) / (speedOfSound[c]*speedOfSound[c]);
                outletMassDensity_[c] = massDensity[c] + massDensityCorrection[c];
//                 }
//             else
//             {
//                 massDensityCorrection[c] = (pressure[c] - outletPressure_) / (speedOfSound[c]*speedOfSound[c]);
//                 outletMassDensity_[c] = massDensity[c] - massDensityCorrection[c];
//             }
            }
            else
            {
                outletMassDensity_[c] = massDensity[c];
            }

            // Liou and Fang, 2000, equation 26 STEP 1
            outletNumberDensity_[c] = outletMassDensity_[c] / molecularMass;

            outletTemperature_[c] = outletPressure_ / (gasConstant*outletMassDensity_[c]);


            outletVelocity_[c] = totalMomentum_[c]/totalMass_[c];

            //velocity correction for each boundary cellI
            if(faceNormalVelocity < VSMALL)
            {
                velocityCorrection[c] = (pressure[c] - outletPressure_) / (massDensity[c]*speedOfSound[c]);
                outletVelocity_[c] += velocityCorrection[c]*n;
            }
            else
            {
                velocityCorrection[c] = (outletPressure_ - pressure[c]) / (massDensity[c]*speedOfSound[c]);
                outletVelocity_[c] -= velocityCorrection[c]*n;
            }

//             Info << "Pressure = " << pressure[c] << endl;
//             Info << "massDensity = " << massDensity[c] << endl;
//             Info << "Temperature = " << translationalTemperature[c] << endl;
        }
    }

//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcLiouFangPressureOutletSpecifiedMolarFraction outlet velocity correction = " << velocityCorrection[(faces_.size()/2)] << endl;
//     }
//
//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcLiouFangPressureOutletSpecifiedMolarFraction outlet velocity = " << outletVelocity_[(faces_.size()/2)] << endl;
//     }

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

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    outletTemperature_[f],
                    mass
                )
            );

                // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (outletVelocity_[f] & -sF/fA )/mostProbableSpeed;

            // From Bird eqn 4.22
            accumulatedParcelsToInsert_[iD][f] +=
                moleFractions_[iD]*
                (
                    fA*outletNumberDensity_[f]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
        }
    }
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::setProperties()
{
    outletPressure_ = readScalar(propsDict_.lookup("outletPressure"));

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcLiouFangPressureOutletSpecifiedMolarFraction::setProperties()")
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
            FatalErrorIn("dsmcLiouFangPressureOutletSpecifiedMolarFraction::dsmcLiouFangPressureOutletSpecifiedMolarFraction()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        moleFractions_[i] = readScalar
        (
            moleFractionsDict.lookup(moleculesReduced[i])
        );
    }

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }

    vibT_.setSize(typeIds_.size());

    forAll(vibT_, m)
    {
        vibT_[m].setSize(nFaces_, 0.0);
    }

    vDof_.setSize(typeIds_.size());

    forAll(vDof_, m)
    {
        vDof_[m].setSize(nFaces_, 0.0);
    }

//     totalVibrationalEnergy_.setSize(typeIds_.size());
//
//     forAll(totalVibrationalEnergy_, m)
//     {
//         totalVibrationalEnergy_[m].setSize(nFaces_, 0.0);
//     }

    nTotalParcelsSpecies_.setSize(typeIds_.size());

    forAll(nTotalParcelsSpecies_, m)
    {
        nTotalParcelsSpecies_[m].setSize(nFaces_, 0.0);
    }
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::setNewBoundaryFields()
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.clear();
    cells_.clear();

    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }

    /***/

    nTimeSteps_ = 0;

    forAll(vibT_, iD)
    {
        vibT_[iD].clear();
        vibT_[iD].setSize(nFaces_, 0.0);
    }
    forAll(accumulatedParcelsToInsert_, iD)
    {
        accumulatedParcelsToInsert_[iD].clear();
        accumulatedParcelsToInsert_[iD].setSize(nFaces_, 0.0);
    }

    forAll(totalVibrationalEnergy_, iD)
    {
        totalVibrationalEnergy_[iD].clear();
        totalVibrationalEnergy_[iD].setSize(nFaces_, 0.0);
    }

    forAll(vDof_, iD)
    {
        vDof_[iD].clear();
        vDof_[iD].setSize(nFaces_, 0.0);
    }

    forAll(nTotalParcelsSpecies_, iD)
    {
        nTotalParcelsSpecies_[iD].clear();
        nTotalParcelsSpecies_[iD].setSize(nFaces_, 0.0);
    }

    // recalculate volume of each boundary cell
    cellVolume_.clear();
    cellVolume_.setSize(nFaces_, 0.0);

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

    // reset the other fields
    UMean_.clear();
    outletVelocity_.clear();
    UCollected_.clear();
    outletNumberDensity_.clear();
    outletMassDensity_.clear();
    outletTemperature_.clear();
    totalMomentum_.clear();
    totalMass_.clear();
    totalMassDensity_.clear();
    totalNumberDensity_.clear();
    totalRotationalEnergy_.clear();
    totalRotationalDof_.clear();
    nTotalParcels_.clear();
    nTotalParcelsInt_.clear();
    mcc_.clear();

    UMean_.setSize(nFaces_, vector::zero);
    outletVelocity_.setSize(nFaces_, vector::zero);
    UCollected_.setSize(nFaces_, vector::zero);
    outletNumberDensity_.setSize(nFaces_, 0.0);
    outletMassDensity_.setSize(nFaces_, 0.0);
    outletTemperature_.setSize(nFaces_, 300.0);
    totalMomentum_.setSize(nFaces_, vector::zero);
    totalMass_.setSize(nFaces_, 0.0);
    totalMassDensity_.setSize(nFaces_, 0.0);
    totalNumberDensity_.setSize(nFaces_, 0.0);
    totalRotationalEnergy_.setSize(nFaces_, 0.0);
    totalRotationalDof_.setSize(nFaces_, 0.0);
    nTotalParcels_.setSize(nFaces_, 0.0);
    nTotalParcelsInt_.setSize(nFaces_, 0.0);
    mcc_.setSize(nFaces_, 0.0);
}



} // End namespace Foam

// ************************************************************************* //
