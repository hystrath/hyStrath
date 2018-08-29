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

The wall temperature is a function of the heat flux from incident parcels.

The radiation equation Q = epsilon*sigma*T^4 is used to calculate the wall temperature.

This class gives a value of temperature for each face on the patch.

\*---------------------------------------------------------------------------*/

#include "dsmcFixedHeatFluxFieldPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcFixedHeatFluxFieldPatch, 0);

addToRunTimeSelectionTable(dsmcPatchBoundary, dsmcFixedHeatFluxFieldPatch, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFixedHeatFluxFieldPatch::dsmcFixedHeatFluxFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    wallTemperature_
    (
        IOobject
        (
            patchName_ + word("_wallTemperature"),
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    EcTot_(nFaces_, 0.0),
    EcTotSum_(nFaces_, 0.0),
    newWallTemperature_(nFaces_, 0.0),
    desiredHeatFlux_(),
    relaxationFactor_(),
    stepCounter_(0),
    nSamples_(0),
    nSamplesAverage_(0),
    referenceCp_(),
    referenceRho_(),
    referenceU_(),
    referenceTemperature_(),
    resetAtOutput_()
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();

    forAll(newWallTemperature_, f)
    {
        newWallTemperature_[f] = temperature_;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFixedHeatFluxFieldPatch::~dsmcFixedHeatFluxFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcFixedHeatFluxFieldPatch::initialConfiguration()
{}


void dsmcFixedHeatFluxFieldPatch::calculateProperties()
{
    stepCounter_++;
    nSamplesAverage_++;

    const scalar deltaT = mesh_.time().deltaTValue(); //TODO cloud_.deltaTValue(p.cell());
    
    forAll(EcTotSum_, f)
    {
        EcTotSum_[f] += EcTot_[f];
        EcTot_[f] = 0.0;
    }
    
    if(stepCounter_ >= nSamples_)
    {
        forAll(EcTotSum_, f)
        { 
            if(fabs(EcTotSum_[f]) > VSMALL) // zero face temperature not allowed!
            {
                const label& faceI = faces_[f];
                scalar fA = mag(mesh_.faceAreas()[faceI]);

                scalar heatFlux = EcTotSum_[f]/(deltaT*nSamplesAverage_*fA);
                
                scalar normalisedDesiredHeatFlux = desiredHeatFlux_ / (referenceCp_*referenceRho_*referenceU_*referenceTemperature_);
                
                scalar normalisedHeatFlux = heatFlux / (referenceCp_*referenceRho_*referenceU_*referenceTemperature_);
                
                scalar oldWallTemperature = newWallTemperature_[f];
                
                scalar deltaWallTemperature = relaxationFactor_*(normalisedHeatFlux - normalisedDesiredHeatFlux)*oldWallTemperature;
                
                newWallTemperature_[f] = oldWallTemperature + deltaWallTemperature;
                    
                if(newWallTemperature_[f] < VSMALL)
                {
                    newWallTemperature_[f] = temperature_;
                }
            }                                    
        }
        
        label wppIndex = patchId_;
        
        forAll(newWallTemperature_, f)
        {
            if(newWallTemperature_[f] > VSMALL)
            {
                wallTemperature_.boundaryFieldRef()[wppIndex][f] = newWallTemperature_[f];
            }
        }
        
        IOdictionary newBoundariesDict
        (
            IOobject
            (
                "boundariesDict",
                time_.system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        
        PtrList<entry> patchBoundaryList(newBoundariesDict.lookup("dsmcPatchBoundaries"));
        
        List< autoPtr<dsmcPatchBoundary> > patchBoundaryModels;
        
        patchBoundaryModels.setSize(patchBoundaryList.size());
        
        forAll(patchBoundaryModels, p)
        {
            const entry& boundaryI = patchBoundaryList[p];
            const dictionary& boundaryIDict = boundaryI.dict();
                        
            const dictionary& patchNameDict = boundaryIDict.subDict("patchBoundaryProperties");
            
            const word& patchName = patchNameDict.lookup("patchName");
          
            if(patchName == patchName_)
            {               
                updateProperties(boundaryIDict);
            }
        }
        
        stepCounter_ = 0.0;
        
        if(resetAtOutput_)
        {
            EcTotSum_ = 0.0;
            nSamplesAverage_ = 0;
        }
    }
}

void dsmcFixedHeatFluxFieldPatch::controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    measurePropertiesBeforeControl(p);

    vector& U = p.U();

    scalar& ERot = p.ERot();
    
    labelList& vibLevel = p.vibLevel();
    
//     label& ELevel = p.ELevel();

    label typeId = p.typeId();
    
    scalar m = cloud_.constProps(typeId).mass();

    scalar preIE = 0.5*m*(U & U) + ERot;
//         + cloud_.constProps(typeId).electronicEnergyList()[ELevel];
        
    forAll(vibLevel, i)
    {
        preIE += vibLevel[i]*physicoChemical::k.value()*cloud_.constProps(typeId).thetaV()[i];
    }
    
    label faceId = findIndex(faces_, p.face());

    vector nw = p.normal();
    nw /= mag(nw);

    // Normal velocity magnitude
    scalar U_dot_nw = U & nw;

    // Wall tangential velocity (flow direction)
    vector Ut = U - U_dot_nw*nw;

    Random& rndGen(cloud_.rndGen());

    while (mag(Ut) < SMALL)
    {
        // If the incident velocity is parallel to the face normal, no
        // tangential direction can be chosen.  Add a perturbation to the
        // incoming velocity and recalculate.

        U = vector
        (
            U.x()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.y()*(0.8 + 0.2*rndGen.sample01<scalar>()),
            U.z()*(0.8 + 0.2*rndGen.sample01<scalar>())
        );

        U_dot_nw = U & nw;

        Ut = U - U_dot_nw*nw;
    }

    // Wall tangential unit vector
    vector tw1 = Ut/mag(Ut);

    // Other tangential unit vector
    vector tw2 = nw^tw1;

    const scalar& T = newWallTemperature_[faceId];
    

    scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof = cloud_.constProps(typeId).rotationalDegreesOfFreedom();
    
    scalar vibrationalDof = cloud_.constProps(typeId).vibrationalDegreesOfFreedom();

    U =
        sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal<scalar>()*tw1
          + rndGen.GaussNormal<scalar>()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
        );

    U += velocity_;

    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);
    
    vibLevel = cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);

    measurePropertiesAfterControl(p, 0.0);
    
    scalar postIE = 0.5*m*(U & U) + ERot;
//         + cloud_.constProps(typeId).electronicEnergyList()[p.ELevel()] ;
        
    forAll(vibLevel, i)
    {
        postIE += vibLevel[i]*physicoChemical::k.value()*cloud_.constProps(typeId).thetaV()[i];
    }
    
//     ELevel = cloud_.equipartitionElectronicLevel
//                     (
//                         T,
//                         cloud_.constProps(typeId).degeneracyList(),
//                         cloud_.constProps(typeId).electronicEnergyList(),
//                         typeId
//                     );
    
    EcTot_[faceId] += cloud_.nParticles(patchId(), faceId)*(preIE - postIE);
}

void dsmcFixedHeatFluxFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{

}

void dsmcFixedHeatFluxFieldPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
    
    setProperties();
}

void dsmcFixedHeatFluxFieldPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("initialTemperature"));
    desiredHeatFlux_ = readScalar(propsDict_.lookup("desiredHeatFlux"));
    relaxationFactor_ = readScalar(propsDict_.lookup("relaxationFactor"));
    nSamples_ = readScalar(propsDict_.lookup("nWallSamples"));
    referenceCp_ = readScalar(propsDict_.lookup("referenceSpecificHeat"));
    referenceRho_ = readScalar(propsDict_.lookup("referenceDensity"));
    referenceU_ = readScalar(propsDict_.lookup("referenceVelocity"));
    referenceTemperature_ = readScalar(propsDict_.lookup("referenceTemperature"));
    resetAtOutput_ = Switch(propsDict_.lookup("resetWallSamples"));
}

} // End namespace Foam

// ************************************************************************* //
