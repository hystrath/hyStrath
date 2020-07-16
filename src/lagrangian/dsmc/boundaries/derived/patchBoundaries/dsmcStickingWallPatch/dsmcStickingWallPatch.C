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

#include "dsmcStickingWallPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcStickingWallPatch, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcStickingWallPatch::setProperties()
{
    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
         FatalErrorIn("dsmcStickingWallPatch::setProperties()")
            << "Cannot have zero typeIds being adsorbed." << nl << "in: "
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

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word& moleculeName(moleculesReduced[i]);

        const label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcStickingWallPatch::setProperties()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    adsorptionProbs_.clear();

    adsorptionProbs_.setSize(typeIds_.size(), 0.0);

    forAll(adsorptionProbs_, i)
    {
        adsorptionProbs_[i] = readScalar
        (
            propsDict_.subDict("adsorptionProbabilities")
                .lookup(moleculesReduced[i])
        );
    }

    residenceTimes_.clear();

    residenceTimes_.setSize(typeIds_.size(), 0.0);

    forAll(residenceTimes_, i)
    {
        if (propsDict_.isDict("residenceTimes"))
        {
            residenceTimes_[i] =
            (
                propsDict_.subDict("residenceTimes")
                    .lookupOrDefault<scalar>(moleculesReduced[i], VGREAT)
            );
        }
        else
        {
            residenceTimes_[i] = VGREAT;
        }
    }

    const scalar saturationLimitPerSquareMeters = propsDict_
        .lookupOrDefault<scalar>("saturationLimit", VGREAT);

    const scalarList& facesArea = mesh_.magSf().boundaryField()[patchId()];

    forAll(saturationLimit_, facei)
    {
        saturationLimit_[facei] =
            saturationLimitPerSquareMeters*facesArea[facei];
    }
}


void dsmcStickingWallPatch::readPatchField()
{
    tmp<volScalarField> tnStuckParticles
    (
        new volScalarField
        (
            IOobject
            (
                "nStuckParticles",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nStuckParticles", dimless, 0.0)
        )
    );

    volScalarField& nStuckParticles = tnStuckParticles.ref();

    nStuckParticles_ = nStuckParticles.boundaryField()[patchId()];

    cloud_.boundaryFluxMeasurements()
        .setBoundarynStuckParticles(patchId(), nStuckParticles_);
}


bool dsmcStickingWallPatch::isNotSaturated(const label facei)
{
    if(saturationLimit(facei) - nStuckParticles(facei) > 1e-6)
    {
        return true;
    }

    return false;
}


void dsmcStickingWallPatch::adsorbParticle
(
    dsmcParcel& p,
    scalar localTemperature
)
{
    //- Particle becomes stuck to wall
    p.setStuck();

    scalarField& wallTemperature = p.stuck().wallTemperature();

    vectorField& wallVectors = p.stuck().wallVectors();

    vector& U = p.U();

    scalar& ERot = p.ERot();

    //- Wall unit normal vector and wall unit tangential vectors
    vector nw, tw1, tw2 = vector::zero;

    dsmcPatchBoundary::calculateWallUnitVectors(p, nw, tw1, tw2);

    const label& typeId = p.typeId();

    const scalar mass = cloud_.constProps(typeId).mass();

    if (localTemperature == 0)
    {
        localTemperature = temperature_;
    }

    // Since the parcel might be stuck to this wall over several time steps, we
    // have to calculate the pre-interaction values with the current RWF as
    // it might be updated in one or several intermediate coordinate system
    // evolve step(s).
    const scalar preIE = p.RWF()*(0.5*mass*(U & U) + ERot
        + cloud_.constProps(typeId).eVib_tot(p.vibLevel())
        + cloud_.constProps(typeId).electronicEnergyList()[p.ELevel()]);

    const vector preIMom = p.RWF()*mass*U;

    wallTemperature[3] = preIE;
    wallVectors[3] = preIMom;

    const scalar& T = localTemperature;

    U = vector::zero;

    const label wppIndex = patchId();
    const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
    const label wppLocalFace = wpp.whichFace(p.face());

    wallTemperature[0] = T;
    wallTemperature[1] = wppIndex;
    wallTemperature[2] = wppLocalFace;

    wallVectors[0] = tw1;
    wallVectors[1] = tw2;
    wallVectors[2] = nw;

    //- Increment the number of stuck particles for this face
    nStuckParticles_[wppLocalFace] += p.RWF() * cloud_.coordSystem().dtModel()
        .nParticles(wppIndex, wppLocalFace);
}


void dsmcStickingWallPatch::testForDesorption(dsmcParcel& p)
{
    const label typeId = p.typeId();
    const label iD = findIndex(typeIds_, typeId);

    if(iD != -1)
    {
        const scalar deltaT = cloud_.deltaTValue(p.cell());

        Random& rndGen = cloud_.rndGen();

        //- Calculate the probability that this particle will be released
        if(rndGen.sample01<scalar>() < deltaT/residenceTimes_[iD])
        {
            const scalar mass = cloud_.constProps(typeId).mass();

            vector& U = p.U();

            U = sqrt
                (
                    physicoChemical::k.value()
                  * p.stuck().wallTemperature()[0]
                  / mass
                )
              * (
                    rndGen.GaussNormal<scalar>()*p.stuck().wallVectors()[0]
                  + rndGen.GaussNormal<scalar>()*p.stuck().wallVectors()[1]
                  - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))
                        *p.stuck().wallVectors()[2]
                );

            //- Decrement the number of stuck particles for this face
            const label& wppIndex = p.stuck().wallTemperature()[1];
            const label& wppLocalFace = p.stuck().wallTemperature()[2];
            nStuckParticles_[wppLocalFace] -= p.RWF() * cloud_.coordSystem()
                .dtModel().nParticles(wppIndex, wppLocalFace);

            measurePropertiesAfterDesorption(p);

            p.deleteStuck();
        }
    }
}


void dsmcStickingWallPatch::measurePropertiesAfterDesorption
(
    dsmcParcel& p,
    const scalar& heatOfReaction
)
{
    if(measurePropertiesAtWall_)
    {
        const label wppIndex = p.stuck().wallTemperature()[1];
        const polyPatch& wpp = mesh_.boundaryMesh()[wppIndex];
        const label wppLocalFace = p.stuck().wallTemperature()[2];

        const scalar fA = mag(wpp.faceAreas()[wppLocalFace]);

        const scalar deltaT = cloud_.deltaTValue(p.cell());

        const dsmcParcel::constantProperties&
            constProps(cloud_.constProps(p.typeId()));

        const scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        const scalar U_dot_nw = p.U() & nw;

        const vector Ut = p.U() - U_dot_nw*nw;

        const scalar rwfDivMagUnfADt =
            p.RWF()/max(mag(U_dot_nw)*fA*deltaT, SMALL);

        cloud_.boundaryFluxMeasurements()
            .rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += rwfDivMagUnfADt;

        if(constProps.rotationalDegreesOfFreedom() > 0)
        {
            cloud_.boundaryFluxMeasurements()
                .rhoNIntBF()[p.typeId()][wppIndex][wppLocalFace] +=
                    rwfDivMagUnfADt;
        }

        if(constProps.nElectronicLevels() > 1)
        {
            cloud_.boundaryFluxMeasurements()
                .rhoNElecBF()[p.typeId()][wppIndex][wppLocalFace] +=
                    rwfDivMagUnfADt;
        }

        cloud_.boundaryFluxMeasurements()
            .rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*rwfDivMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .linearKEBF()[p.typeId()][wppIndex][wppLocalFace] +=
                0.5*m*(p.U() & p.U())*rwfDivMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .mccSpeciesBF()[p.typeId()][wppIndex][wppLocalFace] +=
                m*(p.U() & p.U())*rwfDivMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .momentumBF()[p.typeId()][wppIndex][wppLocalFace] +=
                m*Ut*rwfDivMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .rotationalEBF()[p.typeId()][wppIndex][wppLocalFace] +=
                p.ERot()*rwfDivMagUnfADt;

        cloud_.boundaryFluxMeasurements()
            .rotationalDofBF()[p.typeId()][wppIndex][wppLocalFace] +=
                constProps.rotationalDegreesOfFreedom()*rwfDivMagUnfADt;

        const scalar EVibP_tot = constProps.eVib_tot(p.vibLevel());

        cloud_.boundaryFluxMeasurements()
            .vibrationalEBF()[p.typeId()][wppIndex][wppLocalFace] +=
                EVibP_tot*rwfDivMagUnfADt;

        forAll(p.vibLevel(), mode)
        {
            cloud_.boundaryFluxMeasurements()
                .evmsBF()[p.typeId()][mode][wppIndex][wppLocalFace] +=
                    constProps.eVib_m(mode, p.vibLevel()[mode])
                    * rwfDivMagUnfADt;
        }

        cloud_.boundaryFluxMeasurements()
            .electronicEBF()[p.typeId()][wppIndex][wppLocalFace] +=
                constProps.electronicEnergyList()[p.ELevel()]*rwfDivMagUnfADt;

        // post-interaction energy
        scalar postIE = p.RWF()*(0.5*m*(p.U() & p.U()) + p.ERot() + EVibP_tot
            + constProps.electronicEnergyList()[p.ELevel()]);

        // post-interaction momentum
        const vector postIMom = p.RWF()*m*p.U();

        // Note: these do include the p.RWF() of the parcel at the time it was
        // adsorbed on the wall!
        const scalar preIE = p.stuck().wallTemperature()[3];
        const vector preIMom = p.stuck().wallVectors()[3];

        // Note: Do _not_ use the cloud._nParticles() command here because it
        // assumes the wrong RWF. Because this calculation can happen at
        // completely different time steps the parcel RWF might have been up-
        // dated in the mean time. As the parcel might have been cloned in the
        // mean time (note: the new parcel would then also be adsorbed), we
        // have to use the RWF of the impinging parcel in the pre-interaction
        // calculations and the RWF of the leaving parcel in the post-
        // interaction parcels. The heat of reaction is calculated with the
        // post-interaction RWFs.
        const scalar nParticle = cloud_.coordSystem().dtModel()
            .nParticles(wppIndex, wppLocalFace);

        const scalar deltaQ = nParticle
            *(preIE - postIE
                + (p.RWF()*heatOfReaction*physicoChemical::k.value())
            )
            /(deltaT*fA);

        const vector deltaFD = nParticle*(preIMom - postIMom)/(deltaT*fA);

        cloud_.boundaryFluxMeasurements()
            .qBF()[p.typeId()][wppIndex][wppLocalFace] += deltaQ;

        cloud_.boundaryFluxMeasurements()
            .fDBF()[p.typeId()][wppIndex][wppLocalFace] += deltaFD;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcStickingWallPatch::dsmcStickingWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    adsorptionProbs_(),
    residenceTimes_(),
    nStuckParticles_(mesh_.boundaryMesh()[patchId()].size(), 0.0),
    saturationLimit_(mesh_.boundaryMesh()[patchId()].size(), 0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    readPatchField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcStickingWallPatch::~dsmcStickingWallPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcStickingWallPatch::initialConfiguration()
{}


void dsmcStickingWallPatch::calculateProperties()
{}


void dsmcStickingWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcStickingWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
