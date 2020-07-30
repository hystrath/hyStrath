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

#include "dsmcWangPressureInlet.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcWangPressureInlet, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcWangPressureInlet, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcWangPressureInlet::dsmcWangPressureInlet
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
    inletPressure_(),
    inletTemperature_(),
    n_(),
    cellVolume_(faces_.size(), scalar(0.0)),
    accumulatedParcelsToInsert_(),
    inletVelocity_(faces_.size(), vector::zero),
    totalMomentum_(faces_.size(), vector::zero),
    totalMass_(faces_.size(), scalar(0.0)),
    nTotalParcels_(faces_.size(), scalar(0.0)),
    nTimeSteps_(scalar(0.0)),
    mcc_(faces_.size(), scalar(0.0)),
    UMean_(faces_.size(), vector::zero),
    UCollected_(faces_.size(), vector::zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    // calculate required number density at inlet boundary
    // equation 32, Liou and Fang (2000)
    n_ = inletPressure_ / (physicoChemical::k.value()*inletTemperature_);

    //get volume of each boundary cell
    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcWangPressureInlet::~dsmcWangPressureInlet()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcWangPressureInlet::initialConfiguration()
{}

void dsmcWangPressureInlet::calculateProperties()
{

}

void dsmcWangPressureInlet::controlParcelsBeforeMove()
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
            //rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //Normal unit vector *negative* so normal is pointing into the
            //domain
            vector n = sF;
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            //flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            /* -----------------------------------------------------------------------------*/

            //generate Poisson distributed random number of particles to insert
            //see Tysanner & Garcia, Int. J. Numer. Meth. Fluids 2050; 00:1–12
            //this eliminates non-equilibrium behaviour that does not exist in the corresponding physical system

//             label k = 0;
//
//             const scalar target= exp(-accumulatedParcelsToInsert_[iD][f]);
//
//             scalar p = rndGen.sample01<scalar>();
//
//             while (p > target)
//             {
//                 p *= rndGen.sample01<scalar>();
//                 k += 1;
//             }
//
//             label nParcelsToInsert = k;
//
            /* -----------------------------------------------------------------------------*/

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
}

void dsmcWangPressureInlet::controlParcelsBeforeCollisions()
{

}

void dsmcWangPressureInlet::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    const scalar sqrtPi = sqrt(pi);

    scalar molecularMass = 0.0;
    scalar molarconstantPressureSpecificHeat = 0.0;
    scalar molarconstantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label& typeId_ = typeIds_[iD];

        molecularMass += cloud_.constProps(typeId_).mass()*moleFractions_[iD];
        molarconstantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId_).rotationalDegreesOfFreedom())*moleFractions_[iD];
        molarconstantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId_).rotationalDegreesOfFreedom())*moleFractions_[iD];
    }

     // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma = molarconstantPressureSpecificHeat/molarconstantVolumeSpecificHeat;

    vectorField momentum(faces_.size(), vector::zero);
    vectorField UCollected(faces_.size(), vector::zero);
    scalarField mass(faces_.size(), scalar(0.0));
    scalarField nParcels(faces_.size(), scalar(0.0));
    scalarField mcc(faces_.size(), scalar(0.0));
    scalarField translationalTemperature(faces_.size(), scalar(0.0));
    scalarField numberDensity(faces_.size(), scalar(0.0));
    scalarField massDensity(faces_.size(), scalar(0.0));
    scalarField pressure(faces_.size(), scalar(0.0));
    scalarField speedOfSound(faces_.size(), scalar(0.0));
    scalarField velocityCorrection(faces_.size(), scalar(0.0));

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];

        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[celli];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            label iD = findIndex(typeIds_, p->typeId());

            if(iD != -1)
            {
                momentum[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*p->U();
                mass[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass();
                mcc[c] += cloud_.nParticles(celli)*cloud_.constProps(p->typeId()).mass()*mag(p->U())*mag(p->U());
                nParcels[c] += 1.0;
                UCollected[c] += p->U();
            }
        }

        nTotalParcels_[c] += nParcels[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        mcc_[c] += mcc[c];

        UCollected_[c] += UCollected[c];

        massDensity[c] = totalMass_[c] / (cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if(nTotalParcels_[c] > 1)
        {
            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            translationalTemperature[c] = (1.0/(3.0*physicoChemical::k.value()))
                *(
                    ((mcc_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))))
                    - (
                        (mass[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))
                        )*mag(UMean_[c])*mag(UMean_[c]))
                );

            if(translationalTemperature[c] < VSMALL)
            {
                translationalTemperature[c] = 300.00;
            }

            pressure[c] = numberDensity[c]*physicoChemical::k.value()*translationalTemperature[c];


            speedOfSound[c] = sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n /= mag(n);

//             scalar faceNormalVelocity = (n & inletVelocity_[c]);

            inletVelocity_[c] = totalMomentum_[c]/totalMass_[c];

            if(nTimeSteps_ > 100)
            {
            //velocity correction for each boundary cellI
//             if(faceNormalVelocity < VSMALL)
//             {
                velocityCorrection[c] = (pressure[c] - inletPressure_) / (massDensity[c]*speedOfSound[c]);
                inletVelocity_[c] += velocityCorrection[c]*n;
//             }
//             else
//             {
//                 velocityCorrection[c] = (inletPressure_ - pressure[c]) / (massDensity[c]*speedOfSound[c]);
//                 inletVelocity_[c] -= velocityCorrection[c]*n;
//             }
            }
        }
    }

//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcWangPressureInlet inlet velocity correction = " << velocityCorrection[(faces_.size()/2)] << endl;
//         Pout << "dsmcWangPressureInlet inlet velocity = " << inletVelocity_[(faces_.size()/2)] << endl;
//     }

    // compute number of parcels to insert
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label& typeId = typeIds_[iD];
        scalar mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT = cloud_.deltaTValue(mesh_.boundaryMesh()[patchId_].faceCells()[faceI]);

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
                moleFractions_[iD]*
                (
                    fA*n_*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
        }
    }
}

void dsmcWangPressureInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcWangPressureInlet::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcWangPressureInlet::setProperties()
{
    inletPressure_ = readScalar(propsDict_.lookup("inletPressure"));

    inletTemperature_ = readScalar(propsDict_.lookup("inletTemperature"));

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcWangPressureInlet::setProperties()")
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
            FatalErrorIn("dsmcWangPressureInlet::dsmcWangPressureInlet()")
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

    moleFractions_.clear();  //TODO VINCENT molefractions is not a list of scalarFields!

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
}


void dsmcWangPressureInlet::setNewBoundaryFields()
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

    forAll(accumulatedParcelsToInsert_, iD)
    {
        accumulatedParcelsToInsert_[iD].clear();

        accumulatedParcelsToInsert_[iD].setSize(nFaces_, 0.0);
    }

    // recalculate volume of each boundary cell
    cellVolume_.clear();
    cellVolume_.setSize(nFaces_, 0.0);

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

    // reset other fields
    inletVelocity_.clear();
    totalMomentum_.clear();
    totalMass_.clear();
    nTotalParcels_.clear();
    mcc_.clear();
    UMean_.clear();
    UCollected_.clear();

    inletVelocity_.setSize(nFaces_, vector::zero);
    totalMomentum_.setSize(nFaces_, vector::zero);
    totalMass_.setSize(nFaces_, 0.0);
    nTotalParcels_.setSize(nFaces_, 0.0);
    mcc_.setSize(nFaces_, 0.0);
    UMean_.setSize(nFaces_, vector::zero);
    UCollected_.setSize(nFaces_, vector::zero);
}


} // End namespace Foam

// ************************************************************************* //
