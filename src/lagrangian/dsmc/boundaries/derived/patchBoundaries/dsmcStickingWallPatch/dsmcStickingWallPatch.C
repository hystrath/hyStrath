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
    
    residenceTime_.clear();

    residenceTime_.setSize(typeIds_.size(), 0.0);
    
    forAll(residenceTime_, i)
    {
        if (propsDict_.isDict("residenceTimes"))
        {
            residenceTime_[i] =
            (
                propsDict_.subDict("residenceTimes")
                    .lookupOrDefault<scalar>(moleculesReduced[i], VGREAT)
            );
        }
        else
        {
            residenceTime_[i] = VGREAT;
        }
    }
    
    const scalar saturationLimitPerSquareMeters = propsDict_
        .lookupOrDefault<scalar>("saturationLimit", VGREAT);
        
    const scalarList& facesArea = mesh_.magSf().boundaryField()[patchId()];
    
    forAll(saturationLimit_, facei)
    {
        saturationLimit_[facei] = 
            saturationLimitPerSquareMeters*facesArea[facei]/cloud_.nParticle();
    }
}


bool dsmcStickingWallPatch::isNotSaturated(const label facei)
{
    if(nStuckParcels(facei) < saturationLimit(facei))
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
    
    scalar preIE = 0.0;
    vector preIMom = vector::zero;
    
    const label& typeId = p.typeId();
    
    const scalar mass = cloud_.constProps(typeId).mass();
    
    if (localTemperature == 0)
    {
        localTemperature = temperature_;
    }
    
    preIE = 0.5*mass*(U & U) + ERot + 
        cloud_.constProps(typeId).electronicEnergyList()[p.ELevel()];
    
    forAll(p.vibLevel(), i)
    {
        preIE += p.vibLevel()[i]
            *cloud_.constProps(typeId).thetaV()[i]
            *physicoChemical::k.value();
    }
    
    preIMom = mass*U;
    
    wallTemperature[3] = preIE;
    wallVectors[3] = preIMom;
    
    const scalar& T = localTemperature;
    
    //Random& rndGen = cloud_.rndGen();
    
    U = vector::zero; 
       /*SMALL*sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal()*tw1
          + rndGen.GaussNormal()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))*nw
        );*/
    
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
    nStuckParcels_[wallTemperature[2]]++;
}


void dsmcStickingWallPatch::testForDesorption(dsmcParcel& p)
{ 
    //- Calculate the probability that this particle will be released
    const label& iD = p.typeId(); //findIndex(residenceTime_, p.typeId());
    
    /*Info << "p.typeId() " << p.typeId() << endl;
    Info << "iD " << iD << endl;*/
    
    if(iD != -1)
    {
        const scalar& deltaT = mesh_.time().deltaTValue();

        Random& rndGen = cloud_.rndGen();

        if(rndGen.scalar01() < deltaT/residenceTime_[iD])
        {
            const scalar mass = cloud_.constProps(iD).mass();
            
            vector& U = p.U();
            
            U = sqrt
                (
                    physicoChemical::k.value()
                  * p.stuck().wallTemperature()[0]
                  / mass
                )
              * (
                    rndGen.GaussNormal()*p.stuck().wallVectors()[0]
                  + rndGen.GaussNormal()*p.stuck().wallVectors()[1]
                  - sqrt(-2.0*log(max(1 - rndGen.scalar01(), VSMALL)))
                        *p.stuck().wallVectors()[2]
                );
                
            //- Decrement the number of stuck particles for this face
            nStuckParcels_[p.stuck().wallTemperature()[2]]--;
    
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

        const scalar deltaT = mesh_.time().deltaTValue();

        const dsmcParcel::constantProperties& 
            constProps(cloud_.constProps(p.typeId()));

        const scalar m = constProps.mass();

        vector nw = wpp.faceAreas()[wppLocalFace];
        nw /= mag(nw);

        const scalar U_dot_nw = p.U() & nw;

        const vector Ut = p.U() - U_dot_nw*nw;

        const scalar invMagUnfA = 1/max(mag(U_dot_nw)*fA, SMALL);

        cloud_.boundaryFluxMeasurements()
            .rhoNBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfA;
            
        if(constProps.rotationalDegreesOfFreedom() > 0)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNIntBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfA; 
        }
        
        if(constProps.numberOfElectronicLevels() > 1)
        {
           cloud_.boundaryFluxMeasurements()
              .rhoNElecBF()[p.typeId()][wppIndex][wppLocalFace] += invMagUnfA;
        }

        cloud_.boundaryFluxMeasurements()
            .rhoMBF()[p.typeId()][wppIndex][wppLocalFace] += m*invMagUnfA;
        
        cloud_.boundaryFluxMeasurements()
            .linearKEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                0.5*m*(p.U() & p.U())*invMagUnfA;
        
        cloud_.boundaryFluxMeasurements()
            .mccSpeciesBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*(p.U() & p.U())*invMagUnfA;
        
        cloud_.boundaryFluxMeasurements()
            .momentumBF()[p.typeId()][wppIndex][wppLocalFace] += 
                m*Ut*invMagUnfA;
        
        cloud_.boundaryFluxMeasurements()
            .rotationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                p.ERot()*invMagUnfA;
        
        cloud_.boundaryFluxMeasurements()
            .rotationalDofBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.rotationalDegreesOfFreedom()*invMagUnfA;

        forAll(p.vibLevel(), i)
        {
            cloud_.boundaryFluxMeasurements()
                .vibrationalEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                    p.vibLevel()[i]*constProps.thetaV()[i]
                  * physicoChemical::k.value()
                  * invMagUnfA;
        }

        cloud_.boundaryFluxMeasurements()
            .electronicEBF()[p.typeId()][wppIndex][wppLocalFace] += 
                constProps.electronicEnergyList()[p.ELevel()]*invMagUnfA;

        // post-interaction energy
        scalar postIE = 0.5*m*(p.U() & p.U()) + p.ERot() 
            + constProps.electronicEnergyList()[p.ELevel()];

        forAll(p.vibLevel(), i)
        {
           postIE += p.vibLevel()[i]*constProps.thetaV()[i]
              *physicoChemical::k.value();
        }

        // post-interaction momentum
        const vector postIMom = m*p.U();

        const scalar preIE = p.stuck().wallTemperature()[3];
        const vector preIMom = p.stuck().wallVectors()[3];
        
        scalar nParticle = cloud_.nParticle();
        
        const scalar& RWF = cloud_.getRWF_face(wppLocalFace);
        
        nParticle *= RWF;

        const scalar deltaQ = nParticle
            *(preIE - postIE  + (heatOfReaction*physicoChemical::k.value()))
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
    residenceTime_(),
    nStuckParcels_(mesh_.boundaryMesh()[patchId()].size(), 0.0),
    saturationLimit_(mesh_.boundaryMesh()[patchId()].size(), 0.0)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
    
    setProperties();
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
