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

#include "dsmcZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcZone, 0);

addToRunTimeSelectionTable(dsmcField, dsmcZone, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcZone::dsmcZone
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     fields_(t, mesh, "dummy"),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),
    totalVolume_(0.0),
    typeIds_(),
    timeIndex_(0),
    averagingCounter_(0.0),
    
    mols_(0.0),
    dsmcMols_(0.0),
    molsInt_(0.0),
    molsElec_(0.0),
    mass_(0.0),
    mcc_(0.0),
    mom_(vector::zero),
    UCollected_(vector::zero),
    rotationalEMean_(0.0),
    rotationalDofMean_(0.0),

    muu_(0.0),
    muv_(0.0),
    muw_(0.0),
    mvv_(0.0),
    mvw_(0.0),
    mww_(0.0),
    mccu_(0.0),
    mccv_(0.0),
    mccw_(0.0),

    eu_(0.0),
    ev_(0.0),
    ew_(0.0),
    e_(0.0),
    stepIndex_(0),
    speciesMols_(),
    mccSpecies_(),
    electronicETotal_(),
    nParticlesGroundElectronicState_(),
    nParticlesFirstElectronicState_(),
    vDof_(),
    mfp_(),
    mcr_(),
    
    N_(),
    rhoN_(),
    rhoM_(),
    UMean_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    overallTemperature_(),
    scalarPressure_(),
    pField_(),
    tauField_(),
    qField_(),
    qInternalField_(),
    qTranslationalField_(),
    meanFreePath_(),
    meanCollisionRate_(),
    meanCollisionTime_(),
    meanCollisionTimeTimeStepRatio_(),
    Ma_(),
    vibrationalETotal_(),
    
    outputField_(4, true),
    instantaneous_(false),
    averagingAcrossManyRuns_(false)
    
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcZone::dsmcZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    // standard to reading typeIds ------------ 
    const List<word> molecules (propsDict_.lookup("typeIds"));

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

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcInflowPatch::dsmcInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    //-set the total volume
    const labelList& cells = cellZones[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        totalVolume_ += mesh_.cellVolumes()[cellI];
    }

    if (Pstream::parRun())
    {
        reduce(totalVolume_, sumOp<scalar>());
    }

    // ---------------------------------------------------
    
    
//     instantaneous_ = false;
    
    if (propsDict_.found("instantaneous"))
    {
        instantaneous_ = Switch(propsDict_.lookup("instantaneous"));
    }
    
    // instantaneous
    const scalar& deltaT = time_.mdTimeInterval().deltaT();
    scalar writeInterval = readScalar(t.controlDict().lookup("writeInterval"));
    label nBins = label(writeInterval/deltaT);
    nSteps_ = 1;
    
    if(!instantaneous_) // cumulative
    {
        nBins = 1;
        nSteps_ = label(writeInterval/deltaT);        
    }
    
    N_.setSize(nBins, 0.0);
    rhoN_.setSize(nBins, 0.0);
    rhoM_.setSize(nBins, 0.0);
    UMean_.setSize(nBins, vector::zero);
    UCAM_.setSize(nBins, vector::zero);    
    translationalTemperature_.setSize(nBins, 0.0);
    rotationalTemperature_.setSize(nBins, 0.0);
    vibrationalTemperature_.setSize(nBins, 0.0);
    electronicTemperature_.setSize(nBins, 0.0);
    overallTemperature_.setSize(nBins, 0.0);
    scalarPressure_.setSize(nBins, 0.0);
    pField_.setSize(nBins, tensor::zero);
    tauField_.setSize(nBins, tensor::zero);
    qField_.setSize(nBins, vector::zero);
    qInternalField_.setSize(nBins, vector::zero);
    qTranslationalField_.setSize(nBins, vector::zero);
    vibrationalModeTemperatures_.setSize(nBins, vector::zero);
    meanFreePath_.setSize(nBins, 0.0);
    meanCollisionRate_.setSize(nBins, 0.0);
    meanCollisionTime_.setSize(nBins, 0.0);
    meanCollisionTimeTimeStepRatio_.setSize(nBins, 0.0);
    Ma_.setSize(nBins, 0.0);
    
    speciesMols_.setSize(typeIds_.size(), 0.0);
    mccSpecies_.setSize(typeIds_.size(), 0.0);
    vibrationalETotal_.setSize(typeIds_.size());
    electronicETotal_.setSize(typeIds_.size(), 0.0);
    nParticlesGroundElectronicState_.setSize(typeIds_.size(), 0.0);
    nParticlesFirstElectronicState_.setSize(typeIds_.size(), 0.0);
    vDof_.setSize(typeIds_.size(), 0.0);
    mfp_.setSize(typeIds_.size(), 0.0);
    mcr_.setSize(typeIds_.size(), 0.0);
    
    if (propsDict_.found("averagingAcrossManyRuns"))
    {
        averagingAcrossManyRuns_ = Switch(propsDict_.lookup("averagingAcrossManyRuns"));
        
        // read in stored data from dictionary
        if(averagingAcrossManyRuns_)
        {
            Info << nl << "Averaging across many runs initiated." << nl << endl;

            readIn();
        }         
    }
        
    // choice of measurement property to output

    if (propsDict_.found("outputProperties"))
    {

        const List<word> measurements (propsDict_.lookup("outputProperties"));

        DynamicList<word> propertyNames(0);

        forAll(measurements, i)
        {
            const word& propertyName(measurements[i]);
    
            if(findIndex(propertyNames, propertyName) == -1)
            {
                propertyNames.append(propertyName);
            }
        }

        propertyNames.shrink();

        if(findIndex(propertyNames, "density") == -1)
        {
            outputField_[0] = false;
        }

        if(findIndex(propertyNames, "velocity") == -1)
        {
            outputField_[1] = false;
        }

        if(findIndex(propertyNames, "temperature") == -1)
        {
            outputField_[2] = false;
        }

        if(findIndex(propertyNames, "pressure") == -1)
        {
            outputField_[3] = false;
        }

        // check for user:
        forAll(propertyNames, i)
        {
            const word& propertyName(propertyNames[i]);

            if
            (
                (propertyName != "density") &&
                (propertyName != "velocity") &&
                (propertyName != "temperature") &&
                (propertyName != "pressure")
            )
            {    
                FatalErrorIn("dsmcZone::dsmcZone()")
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << time_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);            
            }
        }
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcZone::~dsmcZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcZone::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "zone_"+fieldName_+"_"+regionName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("mols", mols_);
    dict.readIfPresent("dsmcMols", dsmcMols_);
    dict.readIfPresent("molsInt", molsInt_);
    dict.readIfPresent("molsElec", molsElec_);
    dict.readIfPresent("mass", mass_);
    dict.readIfPresent("mcc", mcc_);
    dict.readIfPresent("mccSpecies", mccSpecies_);
    dict.readIfPresent("mom", mom_);
    dict.readIfPresent("UCollected", UCollected_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
    
    dict.readIfPresent("muu", muu_);    
    dict.readIfPresent("muv", muv_);  
    dict.readIfPresent("muw", muw_);  
    dict.readIfPresent("mvv", mvv_);  
    dict.readIfPresent("mvw", mvw_);  
    dict.readIfPresent("mww", mww_);  

    dict.readIfPresent("mccu", mccu_);  
    dict.readIfPresent("mccv", mccv_);  
    dict.readIfPresent("mccw", mccw_);  
    dict.readIfPresent("eu", eu_);  
    dict.readIfPresent("ev", ev_);      
    dict.readIfPresent("ew", ew_);      
    dict.readIfPresent("e", e_);     
    
    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("electronicETotal", electronicETotal_);
    dict.readIfPresent("nParticlesGroundElectronicState", nParticlesGroundElectronicState_);
    dict.readIfPresent("nParticlesFirstElectronicState", nParticlesFirstElectronicState_);
    dict.readIfPresent("speciesMols", speciesMols_);
    
    dict.readIfPresent("averagingCounter", averagingCounter_);
}

void dsmcZone::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "zone_"+fieldName_+"_"+regionName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("mols", mols_);
        dict.add("dsmcMols", dsmcMols_);
        dict.add("molsInt", molsInt_);
        dict.add("molsElec", molsElec_);
        dict.add("mass", mass_);
        dict.add("mcc", mcc_);
        dict.add("mccSpecies", mccSpecies_);
        dict.add("mom", mom_);
        dict.add("UCollected", UCollected_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        
        dict.add("muu", muu_);    
        dict.add("muv", muv_);  
        dict.add("muw", muw_);  
        dict.add("mvv", mvv_);  
        dict.add("mvw", mvw_);  
        dict.add("mww", mww_);  

        dict.add("mccu", mccu_);  
        dict.add("mccv", mccv_);  
        dict.add("mccw", mccw_);  
        dict.add("eu", eu_);  
        dict.add("ev", ev_);      
        dict.add("ew", ew_);      
        dict.add("e", e_);
        
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("electronicETotal", electronicETotal_);
        dict.add("nParticlesGroundElectronicState", nParticlesGroundElectronicState_);
        dict.add("nParticlesFirstElectronicState", nParticlesFirstElectronicState_);
        dict.add("speciesMols", speciesMols_);
        
        dict.add("averagingCounter", averagingCounter_);
        
        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

void dsmcZone::createField()
{
    forAll(vibrationalETotal_, iD)
    {
        vibrationalETotal_[iD].setSize(cloud_.constProps(typeIds_[iD]).vibrationalDegreesOfFreedom(),0.0);
    }
}


void dsmcZone::calculateField()
{
    averagingCounter_ += 1.0;
      
    const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
    const labelList& cells = mesh_.cellZones()[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];
        
        scalar nParticle = cloud_.nParticle();
                    
        if(cloud_.axisymmetric())
        {
            const point& cC = cloud_.mesh().cellCentres()[cellI];
            scalar radius = cC.y();
            
            scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
            
            nParticle *= RWF;
        }

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];
            label iD = findIndex(typeIds_, p->typeId());
            
            if(iD != -1)
            {
                const dsmcParcel::constantProperties& constProp 
                                = cloud_.constProps(p->typeId());
                                
                const scalar& rotationalDof = 
                    cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom();
                    
                const scalar& mass = constProp.mass()*nParticle;
                const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();

                mols_ += nParticle;
                dsmcMols_ += 1.0;
                mass_ += mass;
                mcc_ += mass*mag(p->U())*mag(p->U());
                UCollected_ += p->U();
                mom_ += mass*p->U();
                rotationalEMean_ += p->ERot();
                rotationalDofMean_ += rotationalDof;
                
                muu_ += mass*sqr(p->U().x());
                muv_ += mass*( (p->U().x()) * (p->U().y()) );
                muw_ += mass*( (p->U().x()) * (p->U().z()) );
                mvv_ += mass*sqr(p->U().y());
                mvw_ += mass*( (p->U().y()) * (p->U().z()) );
                mww_ += mass*sqr(p->U().z());
                mccu_ += mass*mag(p->U())*mag(p->U())*(p->U().x());
                mccv_ += mass*mag(p->U())*mag(p->U())*(p->U().y());
                mccw_ += mass*mag(p->U())*mag(p->U())*(p->U().z());
                
                scalarList EVib(cloud_.constProps(typeIds_[iD]).vibrationalDegreesOfFreedom());
                
                forAll(EVib, i)
                {
                    EVib[i] = 0.0;
                }
                
                forAll(EVib, i)
                {
                    EVib[i] = p->vibLevel()[i]*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV()[i];
                }

                eu_ += nParticle*( p->ERot() + gSum(EVib) )*(p->U().x());
                ev_ += nParticle*( p->ERot() + gSum(EVib) )*(p->U().y());
                ew_ += nParticle*( p->ERot() + gSum(EVib) )*(p->U().z());
                e_ += nParticle*( p->ERot() + gSum(EVib) );
                 
                vibrationalETotal_[iD] += EVib;
                electronicETotal_[iD] += electronicEnergies[p->ELevel()];
                speciesMols_[iD] += 1.0;
                mccSpecies_[iD] += mass*mag(p->U())*mag(p->U());
                
                if(rotationalDof > VSMALL)
                {
                    molsInt_ += 1.0;
                }
                
                if(cloud_.constProps(p->typeId()).numberOfElectronicLevels() > 1)
                {
                    molsElec_ += 1.0;
                }
                
                if(p->ELevel() == 0)
                {
                    nParticlesGroundElectronicState_[iD] += 1.0;
                }
                
                if(p->ELevel() == 1)
                {
                    nParticlesFirstElectronicState_[iD] += 1.0;
                }
            }
        }
    }

    stepIndex_++;
    
    if(stepIndex_ >= nSteps_)
    {
        stepIndex_ = 0;
        
        scalar mass = mass_;
        scalar mols = mols_;
        scalar dsmcMols = dsmcMols_;
        scalar molsInt = molsInt_;
        scalar molsElec = molsElec_;
        scalar mcc = mcc_;
        vector mom = mom_;
        vector UCollected = UCollected_;
        scalar rotationalEMean = rotationalEMean_;
        scalar rotationalDofMean = rotationalDofMean_;

        scalar muu = muu_;
        scalar muv = muv_;
        scalar muw = muw_;
        scalar mvv = mvv_;
        scalar mvw = mvw_;
        scalar mww = mww_;
        scalar mccu = mccu_;
        scalar mccv = mccv_;
        scalar mccw = mccw_;

        scalar eu = eu_;
        scalar ev = ev_;
        scalar ew = ew_;
        scalar e = e_;
        
        List<scalarField> vibrationalETotal = vibrationalETotal_;
        scalarField electronicETotal = electronicETotal_;
        scalarField nParticlesGroundElectronicState = nParticlesGroundElectronicState_;
        scalarField nParticlesFirstElectronicState = nParticlesFirstElectronicState_;
        scalarField speciesMols = speciesMols_;
        scalarField mccSpecies = mccSpecies_;
        
        //- parallel communication

        if(Pstream::parRun())
        {
            reduce(mols, sumOp<scalar>());
            reduce(dsmcMols, sumOp<scalar>());
            reduce(molsInt, sumOp<scalar>());
            reduce(molsElec, sumOp<scalar>());
            reduce(mass, sumOp<scalar>());
            reduce(mcc, sumOp<scalar>());
            reduce(mom, sumOp<vector>());                
            reduce(UCollected, sumOp<vector>());
            reduce(rotationalEMean, sumOp<scalar>());
            reduce(rotationalDofMean, sumOp<scalar>());
            
            reduce(muu, sumOp<scalar>());
            reduce(muv, sumOp<scalar>());
            reduce(muw, sumOp<scalar>());
            reduce(mvv, sumOp<scalar>());
            reduce(mvw, sumOp<scalar>());
            reduce(mww, sumOp<scalar>());
            reduce(mccu, sumOp<scalar>());
            reduce(mccv, sumOp<scalar>());
            reduce(mccw, sumOp<scalar>());

            reduce(eu, sumOp<scalar>());
            reduce(ev, sumOp<scalar>());
            reduce(ew, sumOp<scalar>());
            reduce(e, sumOp<scalar>());
            
            forAll(vibrationalETotal, iD)
            {
                reduce(vibrationalETotal[iD], sumOp<scalarList>());
                reduce(electronicETotal[iD], sumOp<scalar>());
                reduce(speciesMols[iD], sumOp<scalar>());
                reduce(mccSpecies[iD], sumOp<scalar>());
                reduce(nParticlesGroundElectronicState[iD], sumOp<scalar>());
                reduce(nParticlesFirstElectronicState[iD], sumOp<scalar>());
            }
        }
        
        const scalar& volume = totalVolume_;
        label n = timeIndex_;
        
        N_[n] = dsmcMols/averagingCounter_;
        rhoN_[n] = mols/(averagingCounter_*volume);
        rhoM_[n] = mass/(averagingCounter_*volume);
        
        if(dsmcMols > 0.0)
        {
            UMean_[n] = UCollected/dsmcMols;
            
            UCAM_[n] = mom/mass;
            
            translationalTemperature_[n] = (1.0/(3.0*physicoChemical::k.value()))
                                            *(
                                                ((mcc/(mols)))
                                                - (
                                                    (mass/(mols)
                                                    )*mag(UMean_[n])*mag(UMean_[n]))
                                            );
                                            
                
            if(rotationalDofMean > VSMALL)
            {
                rotationalTemperature_[n] = (2.0/physicoChemical::k.value())*(rotationalEMean/rotationalDofMean);
            }
            else
            {
                rotationalTemperature_[n] = 0.0;
            }

            tensor p = tensor::zero;

            p.xx() = rhoN_[n]*( muu/(mols) - ((mass/(mols))*UMean_[n].x()*UMean_[n].x()) );
            p.xy() = rhoN_[n]*( muv/(mols) - ((mass/(mols))*UMean_[n].x()*UMean_[n].y()) );
            p.xz() = rhoN_[n]*( muw/(mols) - ((mass/(mols))*UMean_[n].x()*UMean_[n].z()) );
            p.yx() = p.xy();
            p.yy() = rhoN_[n]*( mvv/(mols) - ((mass/(mols))*UMean_[n].y()*UMean_[n].y()) );
            p.yz() = rhoN_[n]*( mvw/(mols) - ((mass/(mols))*UMean_[n].y()*UMean_[n].z()) );
            p.zx() = p.xz();
            p.zy() = p.yz();
            p.zz() = rhoN_[n]*(mww/(mols) - ((mass/(mols))*UMean_[n].z()*UMean_[n].z()));
            
            pField_[n] = p;

            scalarPressure_[n] = (1.0/3.0)*(p.xx() + p.yy() + p.zz());
                                    
            // make reference 
//             tensorField tau(mols.size(), tensor::zero); 
            tensor tau = tensor::zero;
            
            tau = -p[n];
            tau.xx() += scalarPressure_[n];
            tau.yy() += scalarPressure_[n];
            tau.zz() += scalarPressure_[n];
            tauField_[n] = tau;
            
            vector q = vector::zero;
            
            q.x() = rhoN_[n]*(
                                    0.5*(mccu/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].x()
                                    + eu/(mols)
                                    - (e/(mols))*UMean_[n].x()
                                )
                                    - p.xx()*UMean_[n].x()
                                    - p.xy()*UMean_[n].y()
                                    - p.xz()*UMean_[n].z();
                                    
            //terms involving pressure tensor should not be multiplied 
            //by the number density (see Bird corrigendum)
             
            q.y() = rhoN_[n]*(
                                    0.5*(mccv/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].y()
                                    + ev/(mols)
                                    - (e/(mols))*UMean_[n].y()
                                )
                                    - p.yx()*UMean_[n].x()
                                    - p.yy()*UMean_[n].y()
                                    - p.yz()*UMean_[n].z();
            
            q.z() = rhoN_[n]*(
                                    0.5*(mccw/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].z()
                                    + ew/(mols)
                                    - (e/(mols))*UMean_[n].z()
                                )
                                    - p.zx()*UMean_[n].x()
                                    - p.zy()*UMean_[n].y()
                                    - p.zz()*UMean_[n].z();
            
            qField_[n] = q;
            
            vector qInternal = vector::zero;
            
            qInternal.x() = rhoN_[n]*
                                (
                                    + eu/(mols)
                                    - (e/(mols))*UMean_[n].x()
                                );

            qInternal.y() = rhoN_[n]*
                                (
                                    + ev/(mols)
                                    - (e/(mols))*UMean_[n].y()
                                );
            
            qInternal.z() = rhoN_[n]*
                                (
                                    + ew/(mols)
                                    - (e/(mols))*UMean_[n].z()
                                );
            
            qInternalField_[n] = qInternal;
            
            vector qTranslational = vector::zero;
            
            qTranslational.x() = rhoN_[n]*(
                                    0.5*(mccu/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].x()
                                )
                                    - p.xx()*UMean_[n].x()
                                    - p.xy()*UMean_[n].y()
                                    - p.xz()*UMean_[n].z();            
            
            qTranslational.y() = rhoN_[n]*(
                                    0.5*(mccv/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].y()
                                )
                                    - p.yx()*UMean_[n].x()
                                    - p.yy()*UMean_[n].y()
                                    - p.yz()*UMean_[n].z();
            
            qTranslational.z() = rhoN_[n]*(
                                    0.5*(mccw/(mols))
                                    - 0.5*(mcc/(mols))*UMean_[n].z()
                                )
                                    - p.zx()*UMean_[n].x()
                                    - p.zy()*UMean_[n].y()
                                    - p.zz()*UMean_[n].z();
            
            qTranslationalField_[n] = qTranslational;       
            
            // vibrational temperature
            scalar totalvDof = 0.0;
            scalar totalvDofOverall = 0.0;
            scalar vibT = 0.0;
            scalarList degreesOfFreedomSpecies(typeIds_.size(),0.0);
            scalarList vibTID(vibrationalETotal.size(),0.0);
            
            List<scalarList> degreesOfFreedomMode;
            List<scalarList> vibTMode;
            
            degreesOfFreedomMode.setSize(typeIds_.size());
            vibTMode.setSize(typeIds_.size());
            
            forAll(degreesOfFreedomMode, iD)
            {
                degreesOfFreedomMode[iD].setSize(vibrationalETotal.size(), 0.0);
                vibTMode[iD].setSize(vibrationalETotal.size(), 0.0);
            }
            
            forAll(vibrationalETotal, iD)
            {                
                forAll(vibrationalETotal[iD], m)
                {
                    if(vibrationalETotal[iD][m] > VSMALL && speciesMols[iD] > VSMALL)
                    {        
                        scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV()[m];
                             
                        scalarList vibrationalEMean = (vibrationalETotal[iD]/speciesMols[iD]);
                        
                        scalar iMean = 0.0;
                    
                        iMean = vibrationalEMean[m]/(physicoChemical::k.value()*thetaV);

                        vibTMode[iD][m] = thetaV / log(1.0 + (1.0/iMean));

                        degreesOfFreedomMode[iD][m] = (2.0*thetaV/vibTMode[iD][m]) / (exp(thetaV/vibTMode[iD][m]) - 1.0);
                    }
                }
                    
                forAll(vibrationalETotal[iD], m)
                {
                    degreesOfFreedomSpecies[iD] += degreesOfFreedomMode[iD][m];
                }
                
                forAll(vibrationalETotal[iD], m)
                {
                    if(degreesOfFreedomSpecies[iD] > VSMALL)
                    {
                        vibTID[iD] += vibTMode[iD][m]*degreesOfFreedomMode[iD][m]/degreesOfFreedomSpecies[iD];
                    }
                }
                
                totalvDof += degreesOfFreedomSpecies[iD];
                
                scalar fraction = 0.0;
                   
                if(molsInt > VSMALL)
                {
                    fraction = speciesMols[iD]/molsInt;
                }
                
                scalar fractionOverall = speciesMols[iD]/mols;
                
                if(fraction > SMALL)
                {
                    totalvDofOverall += totalvDof*(fractionOverall/fraction);
                }
                
                vibT += vibTID[iD]*fraction;
            }
            
            vibrationalTemperature_[n] = vibT;
            vibrationalModeTemperatures_[n].x() = vibTMode[0][0];
            vibrationalModeTemperatures_[n].y() = vibTMode[0][1];
            vibrationalModeTemperatures_[n].z() = vibTMode[0][2];

//             // vibrational temperature
//             scalar totalvDof = 0.0;
//             scalar totalvDofOverall = 0.0;
//             scalar vibT = 0.0;
//             
//             forAll(vibrationalETotal, iD)
//             {
//                 if(vibrationalETotal[iD] > VSMALL && speciesMols[iD] > VSMALL)
//                 {        
//                     const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
//                     
//                     scalar vibrationalEMean = (vibrationalETotal[iD]/speciesMols[iD]);
//                     
//                     scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                     
//                     scalar fraction = speciesMols[iD]/molsInt;
//                     
//                     scalar fractionOverall = speciesMols[iD]/dsmcMols;
//                     
//                     scalar vibTID = thetaV / log(1.0 + (1.0/iMean));
//                     
//                     vDof_[iD] = fraction*(2.0*thetaV/vibTID) / (exp(thetaV/vibTID) - 1.0);
//                     
//                     totalvDof += vDof_[iD];
//                     
//                     totalvDofOverall += totalvDof*(fractionOverall/fraction);
//                     
//                     vibT += vibTID*fraction;
//                 }
//             }
//             
//             vibrationalTemperature_[n] = vibT;
   
            // electronic temperature
            scalar totalEDof = 0.0;
            scalar elecT = 0.0;
                
            forAll(speciesMols, iD)
            {
                label nElectronicLevels = cloud_.constProps(typeIds_[iD]).numberOfElectronicLevels();
                
                if(nElectronicLevels > 1 && speciesMols[iD] > VSMALL && molsElec > VSMALL)
                {
                    const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
                    const labelList& degeneracies = cloud_.constProps(typeIds_[iD]).degeneracyList();
                    
//                     scalar speciesTransT = (1.0/(3.0*physicoChemical::k.value()))
//                                             *(
//                                                 ((mccSpecies[iD]/(speciesMols[iD]*cloud_.nParticle())))
//                                                 - (
//                                                     (cloud_.constProps(typeIds_[iD]).mass()/(speciesMols[iD]*cloud_.nParticle())
//                                                     )*mag(UMean_[n])*mag(UMean_[n]))
//                                             );
//                     
//                     scalar fraction = speciesMols[iD]/molsElec;
//                     
//                     if(speciesTransT > VSMALL)
//                     {
//                         scalar sum1 = 0.0;
//                         scalar sum2 = 0.0;
//                         
//                         forAll(electronicEnergies, ii)
//                         {
//                             sum1 += degeneracies[ii]*exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
//                             sum2 += degeneracies[ii]*(electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT))
//                                         *exp(-electronicEnergies[ii]/(physicoChemical::k.value()*speciesTransT));
//                         }
//                         
//                         scalar elecTID = (electronicETotal[iD]/(physicoChemical::k.value()*speciesMols[iD]))*(sum1/sum2);
//                         
//                         elecT += fraction*elecTID;
//                         
//                         scalar eDof = (2.0*(electronicETotal[iD]/speciesMols[iD]))/(physicoChemical::k.value()*speciesTransT);
//                         
//                         totalEDof += fraction*eDof;
//                         
//                         Info << "speciesTransT = " << speciesTransT << endl;
//                     }
                    
                    if(nParticlesGroundElectronicState[iD] > VSMALL && nParticlesFirstElectronicState[iD] > VSMALL && ((nParticlesGroundElectronicState[iD]*degeneracies[1]) != (nParticlesFirstElectronicState[iD]*degeneracies[0])))
                    {
                        scalar fraction = speciesMols[iD]/molsElec;

                        scalar elecTID = (electronicEnergies[1]-electronicEnergies[0])
                            /(physicoChemical::k.value()*
                            log((nParticlesGroundElectronicState[iD]*degeneracies[1])/(nParticlesFirstElectronicState[iD]*degeneracies[0])));
    
                        if(elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }

                        scalar eDof = (2.0*(electronicETotal[iD]/speciesMols[iD]))/(physicoChemical::k.value()*elecTID);
    
                        totalEDof += fraction*eDof;
                    }
                }
            }

            electronicTemperature_[n] = elecT;
            
//             Info << "totalEDof = " << totalEDof << endl;
            
            //overallTemperature
        
            scalar nRotDof = 0.0;
                
            if(dsmcMols > VSMALL)
            {
                nRotDof = rotationalDofMean / dsmcMols;
            }
            
            overallTemperature_[n] = ( 
                                    (3.0*translationalTemperature_[n]) 
                                    + (nRotDof*rotationalTemperature_[n]) 
                                    + (totalvDofOverall*vibrationalTemperature_[n])
                                    + (totalEDof*electronicTemperature_[n])
                                ) /
                                (3.0 + nRotDof + totalvDofOverall + totalEDof);

            forAll(mfp_, iD)
            {
                label qspec = 0;
                
                for (qspec=0; qspec<typeIds_.size(); qspec++)
                {
                    scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());
                    
                    scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());
                    
                    scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();
                    
                    scalar reducedMass = (cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass())
                                            / (cloud_.constProps(typeIds_[iD]).mass()+cloud_.constProps(typeIds_[qspec]).mass());
                    
                    if(speciesMols[qspec] > VSMALL && translationalTemperature_[n] > VSMALL)
                    {
                        scalar nDensQ = (cloud_.nParticle()*speciesMols[qspec])/(volume*averagingCounter_);
                        
                        mfp_[iD] += (pi*dPQ*dPQ*nDensQ*pow(273.0/translationalTemperature_[n],omegaPQ-0.5)*sqrt(1.0+massRatio)); //Bird, eq (4.76)
                        
                        mcr_[iD] += (2.0*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalTemperature_[n]/273.0,1.0-omegaPQ)
                                            *sqrt(2.0*physicoChemical::k.value()*273.0/reducedMass)); // Bird, eq (4.74)
                    }
                }
                
                if(mfp_[iD] > VSMALL)
                {
                    mfp_[iD] = 1.0/mfp_[iD];
                }
            }

            meanFreePath_[n] = 0.0;
            meanCollisionRate_[n] = 0.0;
            meanCollisionTime_[n] = 0.0;
            meanCollisionTimeTimeStepRatio_[n] = 0.0;
                                
            forAll(mfp_, iD)
            {
                if(rhoN_[n] > VSMALL)
                {                    
                    scalar nDensP = (cloud_.nParticle()*speciesMols[iD])/(volume*averagingCounter_);
                    
                    meanFreePath_[n] += mfp_[iD]*nDensP/rhoN_[n]; //Bird, eq (4.77)
                    
                    meanCollisionRate_[n] += mcr_[iD]*nDensP/rhoN_[n]; //Bird, eq (1.38)
                }
                else
                {
                    meanFreePath_[n] = GREAT;
                    meanCollisionRate_[n] = 0.0;
                }
            }
            
            const scalar deltaT = mesh_.time().deltaTValue();

            if(meanCollisionRate_[n] > VSMALL)
            {
                meanCollisionTime_[n] = 1.0/meanCollisionRate_[n];
                meanCollisionTimeTimeStepRatio_[n] = meanCollisionTime_[n]/deltaT;
            }
            else
            {
                meanCollisionTime_[n] = GREAT;
                meanCollisionTimeTimeStepRatio_[n] = GREAT;
            }

            mfp_ = scalar(0.0);
            mcr_ = scalar(0.0);
            
            scalar molecularMass = 0.0;
            scalar molarconstantPressureSpecificHeat = 0.0;
            scalar molarconstantVolumeSpecificHeat = 0.0;
            scalar speedOfSound = 0.0;
            scalar gasConstant = 0.0;
            scalar gamma = 0.0;
            
            forAll(mfp_, iD)  
            {
                const label& typeId = typeIds_[iD];

                if(dsmcMols > VSMALL)
                {
                    molecularMass += cloud_.constProps(typeId).mass()*(speciesMols[iD]/dsmcMols);
                    molarconstantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(speciesMols[iD]/dsmcMols);
                    molarconstantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(speciesMols[iD]/dsmcMols);
                }
            } 

            if(molecularMass > VSMALL)
            {
                gasConstant = physicoChemical::k.value()/molecularMass; // R = k/m
            }

            if(molarconstantVolumeSpecificHeat > VSMALL)
            {
                gamma = molarconstantPressureSpecificHeat/molarconstantVolumeSpecificHeat; // gamma = cP/cV
            }
            
            if(translationalTemperature_[n] > VSMALL && gamma > VSMALL && gasConstant > VSMALL)
            {
                speedOfSound = sqrt(gamma*gasConstant*translationalTemperature_[n]);
            }
            
            if(speedOfSound > VSMALL)
            {
                Ma_[n] = mag(UMean_[n])/speedOfSound;
            }
            else
            {
                Ma_[n] = 0.0;
            }
        }

        if(time_.resetFieldsAtOutput())
        {
            //- reset fields
            averagingCounter_ = 0.0;
            
            mols_ = 0.0;
            dsmcMols_ = 0.0;
            molsInt_ = 0.0;
            molsElec_ = 0.0;
            mass_ = 0.0;
            mcc_ = 0.0;
            mom_ = vector::zero;
            UCollected_ = vector::zero;
            rotationalEMean_ = 0.0;
            rotationalDofMean_ = 0.0;
            
            muu_ = 0.0;
            muv_ = 0.0;
            muw_ = 0.0;
            mvv_ = 0.0;
            mvw_ = 0.0;
            mww_ = 0.0;
            mccu_ = 0.0;
            mccv_ = 0.0;
            mccw_ = 0.0;
            mcc_ = 0.0;
            eu_ = 0.0;
            ev_ = 0.0;
            ew_ = 0.0;
            e_ = 0.0;  
            speciesMols_ = 0.0;
            mccSpecies_ = 0.0;
            
            forAll(electronicETotal_, iD)
            {
                electronicETotal_[iD] = 0.0;
                nParticlesGroundElectronicState_[iD] = 0.0;
                nParticlesFirstElectronicState_[iD] = 0.0;
                
                forAll(vibrationalETotal_[iD], j)
                {
                    vibrationalETotal_[iD][j] = 0.0;
                }
            }
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
        
        timeIndex_++;
    }
}

void dsmcZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        
        timeIndex_ = 0;
        
        if(Pstream::master())
        {
//             fileName timePath(runTime.path()/runTime.timeName()/"uniform");

//             scalarField bins = binModel_->binPositions();
//             vectorField vectorBins = binModel_->bins();
            
//             const scalarField& timeField = time_.averagingTimesInOneWriteInterval();
            scalarField timeField(N_.size(), 0.0);
            const scalar& deltaT = time_.mdTimeInterval().deltaT();
            
            forAll(timeField, t)
            {
                timeField[N_.size()-t-1] = runTime.timeOutputValue() - deltaT*t; 
            }

            
            
            // output densities
            if(outputField_[0])
            {
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_N.xy",
                    timeField,
                    N_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_rhoN.xy",
                    timeField,
                    rhoN_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_rhoM.xy",
                    timeField,
                    rhoM_,
                    true
                );
            }

            // output velocities
            if(outputField_[1])
            {
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_U_SAM.xyz",
                    timeField,
                    UMean_,
                    true
                );


                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_U_CAM.xyz",
                    timeField,
                    UCAM_,
                    true
                );
  
            }

            // output temperature
            if(outputField_[2])
            {
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_translationalTemperature.xy",
                    timeField,
                    translationalTemperature_,
                    true
                );
           
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_rotationalTemperature.xy",
                    timeField,
                    rotationalTemperature_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_vibrationalTemperature.xy",
                    timeField,
                    vibrationalTemperature_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_electronicTemperature.xy",
                    timeField,
                    electronicTemperature_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_overallTemperature.xy",
                    timeField,
                    overallTemperature_,
                    true
                );
            }

            // output pressure
            if(outputField_[3])
            {
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_pressureTensor.xyz",
                    timeField,
                    pField_,
                    true
                );

                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_p.xy",
                    timeField,
                    scalarPressure_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_stressTensor.xyz",
                    timeField,
                    tauField_,
                    true
                );

                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_heatFluxVector.xyz",
                    timeField,
                    qField_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_internalHeatFluxVector.xyz",
                    timeField,
                    qInternalField_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_translationalHeatFluxVector.xyz",
                    timeField,
                    qTranslationalField_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_variableHardSphereMeanFreePath.xyz",
                    timeField,
                    meanFreePath_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_meanCollisionRate.xyz",
                    timeField,
                    meanCollisionRate_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_meanCollisionTime.xyz",
                    timeField,
                    meanCollisionTime_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_meanCollisionTimeTimeStepRatio.xyz",
                    timeField,
                    meanCollisionTimeTimeStepRatio_,
                    true
                );
                
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_Ma.xyz",
                    timeField,
                    Ma_,
                    true
                );
            }
        }
    }
}

void dsmcZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //
