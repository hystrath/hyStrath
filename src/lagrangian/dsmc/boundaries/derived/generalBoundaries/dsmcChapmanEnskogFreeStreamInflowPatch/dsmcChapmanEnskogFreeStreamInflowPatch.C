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

#include "dsmcChapmanEnskogFreeStreamInflowPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcChapmanEnskogFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable
(
    dsmcGeneralBoundary,
    dsmcChapmanEnskogFreeStreamInflowPatch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcChapmanEnskogFreeStreamInflowPatch::dsmcChapmanEnskogFreeStreamInflowPatch
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
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    shearStressMax_(),
    heatFluxMax_(),
    shearStressXComponents_(),
    shearStressYComponents_(),
    shearStressZComponents_(),
    heatFlux_(),
    shearStress_(),
    numberDensities_(),
    accumulatedParcelsToInsert_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcChapmanEnskogFreeStreamInflowPatch::~dsmcChapmanEnskogFreeStreamInflowPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcChapmanEnskogFreeStreamInflowPatch::initialConfiguration()
{}

void dsmcChapmanEnskogFreeStreamInflowPatch::calculateProperties()
{

}

void dsmcChapmanEnskogFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    const scalar sqrtPi = sqrt(pi);

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

           const scalar deltaT =
                cloud_.deltaTValue
                (
                    mesh_.boundaryMesh()[patchId_].faceCells()[f]
                );

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
            accumulatedParcelsToInsert_[i][f] += fA*numberDensities_[i]
                *deltaT*mostProbableSpeed
                *(exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta)))
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
        }
    }

    const vector& faceVelocity = velocity_;
    const scalar& faceTranslationalTemperature = translationalTemperature_;
    const scalar& faceRotationalTemperature = rotationalTemperature_;
    const scalar& faceVibrationalTemperature = vibrationalTemperature_;
    const scalar& faceElectronicTemperature = electronicTemperature_;

    scalar pressure = 0.0;

    forAll(typeIds_, m)
    {
        pressure += numberDensities_[m]*translationalTemperature_*physicoChemical::k.value();
    }

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

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                vector U = vector::zero;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                scalar amplitudeParameter;
                scalar gammaS;

                scalar beta = sqrt(mass/(2.0*translationalTemperature_*physicoChemical::k.value()));

                tensor nShearStress = shearStress_/pressure;
                vector nHeatFlux = heatFlux_*2.0*beta/pressure;

                scalar breakdownParameter = max(
                                        mag(shearStressMax_/pressure),
                                        mag(heatFluxMax_*2.0*beta/pressure)
                                    );

                scalar nTries = 0;

                do
                {
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

                    U =
                        sqrt(physicoChemical::k.value()*faceTranslationalTemperature/mass)
                        *(
                            rndGen.GaussNormal<scalar>()*t1
                            + rndGen.GaussNormal<scalar>()*t2
                        )
//                         + (t1 & faceVelocity)*t1
//                         + (t2 & faceVelocity)*t2
                        + mostProbableSpeed*uNormal*n;

                    amplitudeParameter = 1.0 + 30.0*breakdownParameter;

                    vector Unorm = U*beta;

                    gammaS = 1.0 + (nHeatFlux.x()*Unorm.x() + nHeatFlux.y()*Unorm.y() + nHeatFlux.z()*Unorm.z())*(0.4*magSqr(Unorm) - 1.0)
                            - 2.0*(nShearStress.xy()*Unorm.x()*Unorm.y() + nShearStress.xz()*Unorm.x()*Unorm.z() + nShearStress.yz()*Unorm.y()*Unorm.z())
                            - (nShearStress.xx()*Unorm.x()*Unorm.x() + nShearStress.yy()*Unorm.y()*Unorm.y() + nShearStress.zz()*Unorm.z()*Unorm.z());

//                     Info << "nHeatFlux.x() = " << nHeatFlux.x() << endl;
//                     Info << "Unorm.x() = " << Unorm.x() << endl;
//                     Info << "beta = " << beta << endl;
//                     Info << "magSqr(Unorm) = " << magSqr(Unorm) << endl;
//                     Info << "amplitudeParameter = " << amplitudeParameter << endl;
//                     Info << "gammaS = " << gammaS << endl;

//                     nTries += 1;
//                     Info << "nTries = " << nTries << endl;

                }
                while(amplitudeParameter*cloud_.rndGen().sample01<scalar>() >= gammaS);

                U += ((t1 & faceVelocity)*t1 + (t2 & faceVelocity)*t2);

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).nVibrationalModes(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceElectronicTemperature,
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

void dsmcChapmanEnskogFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{

}

void dsmcChapmanEnskogFreeStreamInflowPatch::controlParcelsAfterCollisions()
{
}

void dsmcChapmanEnskogFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcChapmanEnskogFreeStreamInflowPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

//     setProperties();
}



void dsmcChapmanEnskogFreeStreamInflowPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    translationalTemperature_ = readScalar(propsDict_.lookup("translationalTemperature"));
    rotationalTemperature_ = readScalar(propsDict_.lookup("rotationalTemperature"));
    vibrationalTemperature_ = readScalar(propsDict_.lookup("vibrationalTemperature"));

    shearStressXComponents_ = propsDict_.lookup("shearStressXComponents");
    shearStressYComponents_ = propsDict_.lookup("shearStressYComponents");
    shearStressZComponents_ = propsDict_.lookup("shearStressZComponents");
    heatFlux_ = propsDict_.lookup("heatFlux");

    shearStress_.xx() = shearStressXComponents_.x();
    shearStress_.xy() = shearStressXComponents_.y();
    shearStress_.xz() = shearStressXComponents_.z();

    shearStress_.yx() = shearStressYComponents_.x();
    shearStress_.yy() = shearStressYComponents_.y();
    shearStress_.yz() = shearStressYComponents_.z();

    shearStress_.zx() = shearStressZComponents_.x();
    shearStress_.zy() = shearStressZComponents_.y();
    shearStress_.zz() = shearStressZComponents_.z();

    heatFluxMax_ = abs(heatFlux_.x());
    if(abs(heatFlux_.y()) > heatFluxMax_)
    {
        heatFluxMax_ = abs(heatFlux_.y());
    }
    if(abs(heatFlux_.z()) > heatFluxMax_)
    {
        heatFluxMax_ = abs(heatFlux_.z());
    }

    shearStressMax_ = abs(shearStress_.xx());
    if(abs(shearStress_.xy()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.xy());
    }
    if(abs(shearStress_.xz()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.xz());
    }
    if(abs(shearStress_.yx()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.yx());
    }
    if(abs(shearStress_.yy()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.yy());
    }
    if(abs(shearStress_.yz()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.yz());
    }
    if(abs(shearStress_.zx()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.zx());
    }
    if(abs(shearStress_.zy()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.zy());
    }
    if(abs(shearStress_.zz()) > shearStressMax_)
    {
        shearStressMax_ = abs(shearStress_.zz());
    }

    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcChapmanEnskogFreeStreamInflowPatch::dsmcChapmanEnskogFreeStreamInflowPatch()")
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
            FatalErrorIn("dsmcChapmanEnskogFreeStreamInflowPatch::dsmcChapmanEnskogFreeStreamInflowPatch()")
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


void dsmcChapmanEnskogFreeStreamInflowPatch::setNewBoundaryFields()
{
    Info << "dsmcFreeStreamInflowPatch::setNewBoundaryFields()" << endl;

   forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}


} // End namespace Foam

// ************************************************************************* //
