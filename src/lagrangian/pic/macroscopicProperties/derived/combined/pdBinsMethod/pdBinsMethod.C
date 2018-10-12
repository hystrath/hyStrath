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

#include "pdBinsMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(pdBinsMethod, 0);

addToRunTimeSelectionTable(pdField, pdBinsMethod, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pdBinsMethod::pdBinsMethod
(
    Time& t,
    const polyMesh& mesh,
    pdCloud& cloud,
    const dictionary& dict
)
:
    pdField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),
    typeIds_(),

    averagingCounter_(0.0),
    
    volumeOfCellsInBin_(),
    mols_(),
    molsInt_(),
    mass_(),
    mcc_(),
    UCollected_(),
    rotationalEMean_(),
    rotationalDofMean_(),

    muu_(),
    muv_(),
    muw_(),
    mvv_(),
    mvw_(),
    mww_(),
    scalarPressure_(),
    mccu_(),
    mccv_(),
    mccw_(),

    eu_(),
    ev_(),
    ew_(),
    e_(),
    
    speciesMols_(),
    vibrationalETotal_(),
    vDof_(),
    mfp_(),

    binVolume_(),
    N_(),
    M_(),
    rhoN_(),
    rhoM_(),
    averageMass_(),
    UMean_(),
    UCAM_(),
    scalarUMean_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    overallTemperature_(),
    
    pField_(),
    tauField_(),
    qField_(),
    qInternalField_(),
    qTranslationalField_(),
    meanFreePath_(),
    
    outputField_(4, true),
    averagingAcrossManyRuns_(false),
    permeabilityMeasurements_(false),
    rotationalMeasurements_(false),
    wallVelocity_(0.0),
    wallRadius_(0.0) 
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("pdBinsMethod::pdBinsMethod()")
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
            FatalErrorIn("pdBinsMethod::pdBinsMethod()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    
    // ---------------------------------------------------

    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );
    
    const label& nBins = binModel_->nBins();

    volumeOfCellsInBin_.setSize(nBins, 0.0);
    mols_.setSize(nBins, 0.0);
    molsInt_.setSize(nBins, 0.0);
    mass_.setSize(nBins, 0.0);
    mcc_.setSize(nBins, 0.0);
    mom_.setSize(nBins, vector::zero);
    UCollected_.setSize(nBins, vector::zero);
    rotationalEMean_.setSize(nBins, 0.0);
    rotationalDofMean_.setSize(nBins, 0.0);  

    muu_.setSize(nBins, 0.0);
    muv_.setSize(nBins, 0.0);
    muw_.setSize(nBins, 0.0);
    mvv_.setSize(nBins, 0.0);
    mvw_.setSize(nBins, 0.0);
    mww_.setSize(nBins, 0.0);
  
    mccu_.setSize(nBins, 0.0);
    mccv_.setSize(nBins, 0.0);
    mccw_.setSize(nBins, 0.0);

    eu_.setSize(nBins, 0.0);
    ev_.setSize(nBins, 0.0);
    ew_.setSize(nBins, 0.0);
    e_.setSize(nBins, 0.0);
    
    vibrationalETotal_.setSize(nBins);
    speciesMols_.setSize(nBins);
    vDof_.setSize(nBins);
    mfp_.setSize(nBins);
    
    binVolume_.setSize(nBins, 0.0);
    N_.setSize(nBins, 0.0);
    M_.setSize(nBins, 0.0);
    rhoN_.setSize(nBins, 0.0);
    rhoM_.setSize(nBins, 0.0);
    averageMass_.setSize(nBins, 0.0);
    UMean_.setSize(nBins, vector::zero);
    UCAM_.setSize(nBins, vector::zero);  
    scalarUMean_.setSize(nBins, 0.0);
    translationalTemperature_.setSize(nBins, 0.0);
    rotationalTemperature_.setSize(nBins, 0.0);
    vibrationalTemperature_.setSize(nBins, 0.0);
    overallTemperature_.setSize(nBins, 0.0);
    
    scalarPressure_.setSize(nBins, 0.0);
    pField_.setSize(nBins, tensor::zero);
    tauField_.setSize(nBins, tensor::zero);
    qField_.setSize(nBins, vector::zero);
    qInternalField_.setSize(nBins, vector::zero);
    qTranslationalField_.setSize(nBins, vector::zero);
    meanFreePath_.setSize(nBins, 0.0);
    Ma_.setSize(nBins, 0.0);
    
    // outer list is the bins, inner list is type ids
    forAll(vibrationalETotal_, b)
    {
        vibrationalETotal_[b].setSize(typeIds_.size(), 0.0);    
        speciesMols_[b].setSize(typeIds_.size(), 0.0);
        vDof_[b].setSize(typeIds_.size(), 0.0);
        mfp_[b].setSize(typeIds_.size(), 0.0);
    } 
    
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
    
    if (propsDict_.found("permeabilityMeasurements"))
    {
        permeabilityMeasurements_ = Switch(propsDict_.lookup("permeabilityMeasurements"));
        
        // read in stored data from dictionary
        if(permeabilityMeasurements_)
        {
            Info << nl << "Measuring volume of cells in bins; ensure your bins are aligned with the mesh!" << nl << endl;
        }         
    }
    
    if (propsDict_.found("rotationalMeasurements"))
    {
        rotationalMeasurements_ = Switch(propsDict_.lookup("rotationalMeasurements"));
        
        // read in stored data from dictionary
        if(rotationalMeasurements_)
        {
            wallVelocity_ = readScalar(propsDict_.lookup("wallVelocity"));
            
            wallRadius_ = readScalar(propsDict_.lookup("wallRadius"));
            
            if(wallRadius_ < VSMALL)
            {
                FatalErrorIn("pdBinsMethod::pdBinsMethod()")
                << "Wall radius is: " << wallRadius_ << ", must be greater than zero."
                << exit(FatalError);
            }
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
                FatalErrorIn("pdBinsMethod::pdBinsMethod()")
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << time_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);            
            }
        }
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pdBinsMethod::~pdBinsMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pdBinsMethod::readIn()
{
    IOdictionary dict
    (
        IOobject
        (
            "binsMethod_"+fieldName_+"_"+regionName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );   
    
    dict.readIfPresent("volumeOfCellsInBin", volumeOfCellsInBin_);
    dict.readIfPresent("mols", mols_);
    dict.readIfPresent("molsInt", molsInt_);
    dict.readIfPresent("mass", mass_);
    dict.readIfPresent("mcc", mcc_);
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
    dict.readIfPresent("speciesMols", speciesMols_);
    
    dict.readIfPresent("averagingCounter", averagingCounter_);
}

void pdBinsMethod::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "binsMethod_"+fieldName_+"_"+regionName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        
        dict.add("volumeOfCellsInBin", volumeOfCellsInBin_);
        dict.add("mols", mols_);
        dict.add("molsInt", molsInt_);
        dict.add("mass", mass_);
        dict.add("mcc", mcc_);
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
        dict.add("speciesMols", speciesMols_);
        
        dict.add("averagingCounter", averagingCounter_);
        
        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

void pdBinsMethod::createField()
{  
}


void pdBinsMethod::calculateField()
{
    averagingCounter_ += 1.0;
      
    const List< DynamicList<pdParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
    const labelList& cells = mesh_.cellZones()[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        if(averagingCounter_ < SMALL+1.0)
        {
            const vector& cC = mesh_.cellCentres()[cellI];
            
            label bin = binModel_->isPointWithinBin(cC, cellI);
            
            if(bin != -1)
            {
                const scalar& cellVolume = mesh_.cellVolumes()[cellI];
                
                volumeOfCellsInBin_[bin] += cellVolume;
            }
        }
        
        const List<pdParcel*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            pdParcel* p = molsInCell[mIC];
            
            label iD = findIndex(typeIds_, p->typeId());

            const vector& rI = p->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(iD != -1)
                {
                    const pdParcel::constantProperties& constProp 
                                    = cloud_.constProps(p->typeId());
                                    
                    const scalar& rotationalDof = 
                        cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom();
                        
                    const scalar& mass = constProp.mass()*cloud_.nParticle();

                    mols_[n] += 1.0;
                    mass_[n] += mass;
                    mcc_[n] += mass*mag(p->U())*mag(p->U());                   
                    UCollected_[n] += p->U();
                    mom_[n] += mass*p->U();
                    rotationalEMean_[n] += p->ERot();
                    rotationalDofMean_[n] += rotationalDof;
                    
                    muu_[n] += mass*sqr(p->U().x());
                    muv_[n] += mass*( (p->U().x()) * (p->U().y()) );
                    muw_[n] += mass*( (p->U().x()) * (p->U().z()) );
                    mvv_[n] += mass*sqr(p->U().y());
                    mvw_[n] += mass*( (p->U().y()) * (p->U().z()) );
                    mww_[n] += mass*sqr(p->U().z());
                    mccu_[n] += mass*mag(p->U())*mag(p->U())*(p->U().x());
                    mccv_[n] += mass*mag(p->U())*mag(p->U())*(p->U().y());
                    mccw_[n] += mass*mag(p->U())*mag(p->U())*(p->U().z());

                    eu_[n] += cloud_.nParticle()*( p->ERot() + p->EVib() )*(p->U().x());
                    ev_[n] += cloud_.nParticle()*( p->ERot() + p->EVib() )*(p->U().y());
                    ew_[n] += cloud_.nParticle()*( p->ERot() + p->EVib() )*(p->U().z());
                    e_[n] += cloud_.nParticle()*( p->ERot() + p->EVib() );

                    vibrationalETotal_[n][iD] += p->EVib();
                    speciesMols_[n][iD] += 1.0;
                    
                    if(rotationalDof > VSMALL)
                    {
                        molsInt_[n] += 1.0;
                    }
                }
            }
        }
    }
    
    
 
    if(time_.averagingTime())
    {
        scalarField volumeOfCellsInBin = volumeOfCellsInBin_;
        scalarField mass = mass_;
        scalarField molsInt = molsInt_;
        scalarField mols = mols_;
        scalarField mcc = mcc_;
        vectorField mom = mom_;
        vectorField UCollected = UCollected_;
        scalarField rotationalEMean = rotationalEMean_;
        scalarField rotationalDofMean = rotationalDofMean_;

        scalarField muu = muu_;
        scalarField muv = muv_;
        scalarField muw = muw_;
        scalarField mvv = mvv_;
        scalarField mvw = mvw_;
        scalarField mww = mww_;
        scalarField mccu = mccu_;
        scalarField mccv = mccv_;
        scalarField mccw = mccw_;

        scalarField eu = eu_;
        scalarField ev = ev_;
        scalarField ew = ew_;
        scalarField e = e_;
        
        List<scalarField> vibrationalETotal = vibrationalETotal_;
        List<scalarField> speciesMols = speciesMols_;
        
        //- parallel communication

        if(Pstream::parRun())
        {
            forAll(mols, n)
            {
                reduce(volumeOfCellsInBin[n], sumOp<scalar>());
                reduce(mols[n], sumOp<scalar>());
                reduce(molsInt[n], sumOp<scalar>());
                reduce(mass[n], sumOp<scalar>());
                reduce(mcc[n], sumOp<scalar>());
                reduce(mom[n], sumOp<vector>());                
                reduce(UCollected[n], sumOp<vector>());
                reduce(rotationalEMean[n], sumOp<scalar>());
                reduce(rotationalDofMean[n], sumOp<scalar>());
                
                reduce(muu[n], sumOp<scalar>());
                reduce(muv[n], sumOp<scalar>());
                reduce(muw[n], sumOp<scalar>());
                reduce(mvv[n], sumOp<scalar>());
                reduce(mvw[n], sumOp<scalar>());
                reduce(mww[n], sumOp<scalar>());
                reduce(mccu[n], sumOp<scalar>());
                reduce(mccv[n], sumOp<scalar>());
                reduce(mccw[n], sumOp<scalar>());

                reduce(eu[n], sumOp<scalar>());
                reduce(ev[n], sumOp<scalar>());
                reduce(ew[n], sumOp<scalar>());
                reduce(e[n], sumOp<scalar>());
            }
            
            forAll(vibrationalETotal, b)
            {
                forAll(vibrationalETotal[b], n)
                {
                    reduce(vibrationalETotal[b][n], sumOp<scalar>());
                    reduce(speciesMols[b][n], sumOp<scalar>());
                }
            }
        }

        forAll(mols, n)
        {
            scalar volume = binModel_->binVolume(n);
            
            binVolume_[n] = volumeOfCellsInBin[n];
            N_[n] = mols[n]/averagingCounter_;
            M_[n] = mass[n]/averagingCounter_;
            rhoN_[n] = (mols[n]*cloud_.nParticle())/(averagingCounter_*volume);
            rhoM_[n] = mass[n]/(averagingCounter_*volume);
            
            if(mols[n] > 0.0)
            {  
                averageMass_[n] = M_[n]/(N_[n]*cloud_.nParticle());
                
                UMean_[n] = UCollected[n]/mols[n];
                
                if(rotationalMeasurements_)
                {
                    scalarUMean_[n] = wallVelocity_*binModel_->binPositions()[n]/wallRadius_;
                }
                
                UCAM_[n] = mom[n]/mass[n];
                
                translationalTemperature_[n] = (1.0/(3.0*physicoChemical::k.value()))
                                                *(
                                                    ((mcc[n]/(mols[n]*cloud_.nParticle())))
                                                    - (
                                                        (mass[n]/(mols[n]*cloud_.nParticle())
                                                      )*mag(UMean_[n])*mag(UMean_[n]))
                                                );
                                                
                if(rotationalMeasurements_)
                {
                    translationalTemperature_[n] = (1.0/(3.0*physicoChemical::k.value()))
                                                    *(
                                                        ((mcc[n]/(mols[n]*cloud_.nParticle())))
                                                        - (
                                                            (mass[n]/(mols[n]*cloud_.nParticle())
                                                        )*sqr(scalarUMean_[n]))
                                                    );   
                }
                                                
                if(rotationalDofMean[n] > VSMALL)
                {
                    rotationalTemperature_[n] = (2.0/physicoChemical::k.value())*(rotationalEMean[n]/rotationalDofMean[n]);
                }
                else
                {
                    rotationalTemperature_[n] = 0.0;
                }

                tensorField p(mols.size(), tensor::zero);

                p[n].xx() = rhoN_[n]*( muu[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()*UMean_[n].x()) );
                p[n].xy() = rhoN_[n]*( muv[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()*UMean_[n].y()) );
                p[n].xz() = rhoN_[n]*( muw[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()*UMean_[n].z()) );
                p[n].yx() = p[n].xy();
                p[n].yy() = rhoN_[n]*( mvv[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()*UMean_[n].y()) );
                p[n].yz() = rhoN_[n]*( mvw[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()*UMean_[n].z()) );
                p[n].zx() = p[n].xz();
                p[n].zy() = p[n].yz();
                p[n].zz() = rhoN_[n]*(mww[n]/(mols[n]*cloud_.nParticle()) - ((mass[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].z()*UMean_[n].z()));
                
                pField_[n] = p[n];

                scalarPressure_[n] = (1.0/3.0)*(p[n].xx() + p[n].yy() + p[n].zz());
                                        
                // make reference 
                tensorField tau(mols.size(), tensor::zero); 
                
                tau[n] = -p[n];
                tau[n].xx() += scalarPressure_[n];
                tau[n].yy() += scalarPressure_[n];
                tau[n].zz() += scalarPressure_[n];
                tauField_[n] = tau[n];
                
                vectorField q(mols_.size(), vector::zero);
                
                q[n].x() = rhoN_[n]*(
                                        0.5*(mccu[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()
                                        + eu[n]/(mols[n]*cloud_.nParticle())
                                        - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()
                                  )
                                        - p[n].xx()*UMean_[n].x()
                                        - p[n].xy()*UMean_[n].y()
                                        - p[n].xz()*UMean_[n].z();
                                        
                 //terms involving pressure tensor should not be multiplied by the number density (see Bird corrigendum)
                
                
                q[n].y() = rhoN_[n]*(
                                        0.5*(mccv[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()
                                        + ev[n]/(mols[n]*cloud_.nParticle())
                                        - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()
                                  )
                                        - p[n].yx()*UMean_[n].x()
                                        - p[n].yy()*UMean_[n].y()
                                        - p[n].yz()*UMean_[n].z();
                
                q[n].z() = rhoN_[n]*(
                                        0.5*(mccw[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].z()
                                        + ew[n]/(mols[n]*cloud_.nParticle())
                                        - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].z()
                                  )
                                        - p[n].zx()*UMean_[n].x()
                                        - p[n].zy()*UMean_[n].y()
                                        - p[n].zz()*UMean_[n].z();
                
                qField_[n] = q[n];
                
                vectorField qInternal(mols_.size(), vector::zero);
                
                qInternal[n].x() = rhoN_[n]*(
                                            + eu[n]/(mols[n]*cloud_.nParticle())
                                            - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()
                                        );
                
                qInternal[n].y() = rhoN_[n]*(
                                            + ev[n]/(mols[n]*cloud_.nParticle())
                                            - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()
                                        );
                
                qInternal[n].z() = rhoN_[n]*(
                                            + ew[n]/(mols[n]*cloud_.nParticle())
                                            - (e[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].z()
                                        );
                
                qInternalField_[n] = qInternal[n];
                
                
                vectorField qTranslational(mols_.size(), vector::zero);
                
                qTranslational[n].x() = rhoN_[n]*(
                                        0.5*(mccu[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].x()
                                  )
                                        - p[n].xx()*UMean_[n].x()
                                        - p[n].xy()*UMean_[n].y()
                                        - p[n].xz()*UMean_[n].z();
                                        
                qTranslational[n].y() = rhoN_[n]*(
                                        0.5*(mccv[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].y()
                                  )
                                        - p[n].yx()*UMean_[n].x()
                                        - p[n].yy()*UMean_[n].y()
                                        - p[n].yz()*UMean_[n].z();
                
                qTranslational[n].z() = rhoN_[n]*(
                                        0.5*(mccw[n]/(mols[n]*cloud_.nParticle()))
                                        - 0.5*(mcc[n]/(mols[n]*cloud_.nParticle()))*UMean_[n].z()
                                  )
                                        - p[n].zx()*UMean_[n].x()
                                        - p[n].zy()*UMean_[n].y()
                                        - p[n].zz()*UMean_[n].z();
                
                qTranslationalField_[n] = qTranslational[n];
                                        
                const label& nBins = binModel_->nBins();
                
                // vibrational temperature
                scalarField totalvDof(nBins, 0.0);
                scalarField vibT(nBins, 0.0);
                
                forAll(vibrationalETotal[n], iD)
                {
                    if(vibrationalETotal[n][iD] > VSMALL && speciesMols[n][iD] > VSMALL)
                    {        
                        const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
                        
                        scalar vibrationalEMean = (vibrationalETotal[n][iD]/speciesMols[n][iD]);
                        
                        scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
                        
                        scalar fraction = speciesMols[n][iD]/molsInt[n];
                        
                        scalar vibTID = thetaV / log(1.0 + (1.0/iMean));
                        
                        scalar vDof = fraction*(2.0*thetaV/vibTID) / (exp(thetaV/vibTID) - 1.0);
                        
                        totalvDof[n] += vDof;
                        
                        vibT[n] += vibTID*fraction;
                    }
                }
                
                vibrationalTemperature_[n] = vibT[n];
                
                //overallTemperature
        
                scalar nRotDof = 0.0;
                    
                if(mols[n] > VSMALL)
                {
                    nRotDof = rotationalDofMean[n] / mols[n];
                }
                
                overallTemperature_[n] = ( 
                                        (3.0*translationalTemperature_[n]) 
                                        + (nRotDof*rotationalTemperature_[n]) 
                                        + (totalvDof[n]*vibrationalTemperature_[n])
                                    ) /
                                    (3.0 + nRotDof + totalvDof[n]);
                                    
                forAll(mfp_[n], iD)
                {
                    label qspec = 0;
                    
                    for (qspec=0; qspec<typeIds_.size(); qspec++)
                    {
                        scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());
                        scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());
                        scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();
                        
                        if(speciesMols[n][qspec] > VSMALL && translationalTemperature_[n] > VSMALL)
                        {
                            scalar nDensQ = (cloud_.nParticle()*speciesMols[n][qspec])/(volume*averagingCounter_);
                            
                            mfp_[n][iD] += (pi*dPQ*dPQ*nDensQ*pow(273.0/translationalTemperature_[n],omegaPQ-0.5)*sqrt(1.0+massRatio)); //Bird, eq (4.76)
                        }
                    }
                    
                    if(mfp_[n][iD] > VSMALL)
                    {
                        mfp_[n][iD] = 1.0/mfp_[n][iD];
                    }
                }
                
                meanFreePath_[n] = 0.0;
                               
                forAll(mfp_[n], iD)
                {
                    if(rhoN_[n] > VSMALL)
                    {                        
                        scalar nDensP = (cloud_.nParticle()*speciesMols[n][iD])/(volume*averagingCounter_);
                        
                        meanFreePath_[n] += mfp_[n][iD]*nDensP/rhoN_[n]; //Bird, eq (4.77)
                    }
                    else
                    {
                        meanFreePath_[n] = GREAT;
                    }
                }
                
                mfp_[n] = scalar(0.0);
                
                scalar molecularMass = 0.0;
                scalar molarconstantPressureSpecificHeat = 0.0;
                scalar molarconstantVolumeSpecificHeat = 0.0;
                scalar gasConstant = 0.0;
                scalar gamma = 0.0;
                scalar speedOfSound = 0.0;
                
                forAll(mfp_[n], iD)  
                {
                    const label& typeId = typeIds_[iD];

                    if(rhoN_[n] > VSMALL)
                    {
                        molecularMass += cloud_.constProps(typeId).mass()*(speciesMols[n][iD]/mols[n]);
                        molarconstantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(speciesMols[n][iD]/mols[n]);
                        molarconstantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*(speciesMols[n][iD]/mols[n]);
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
        }

        if(time_.resetFieldsAtOutput())
        {
            //- reset fields
            averagingCounter_ = 0.0;
            
            volumeOfCellsInBin_ = 0.0;
            mols_ = 0.0;
            molsInt_ = 0.0;
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
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
    }
}

void pdBinsMethod::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {        
        if(Pstream::master())
        {
            scalarField bins = binModel_->binPositions();
            vectorField vectorBins = binModel_->bins();

            // output densities
            if(outputField_[0])
            {
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_N.xy",
                    bins,
                    N_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_N_3D_pos.xy",
                    vectorBins,
                    N_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoN.xy",
                    bins,
                    rhoN_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoN_3D_pos.xy",
                    vectorBins,
                    rhoN_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoM.xy",
                    bins,
                    rhoM_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoM_3D_pos.xy",
                    vectorBins,
                    rhoM_
                );
            }

            // output velocities
            if(outputField_[1])
            {
                
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_SAM.xyz",
                    bins,
                    UMean_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_SAM_3D_pos.xyz",
                    vectorBins,
                    UMean_
                );    

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_CAM.xyz",
                    bins,
                    UCAM_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_CAM_3D_pos.xyz",
                    vectorBins,
                    UCAM_
                );   
            }

            // output temperature
            if(outputField_[2])
            {
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_translationalTemperature.xy",
                    bins,
                    translationalTemperature_
                );
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_translationalTemperature_3D_pos.xyz",
                    vectorBins,
                    translationalTemperature_
                );                
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rotationalTemperature.xy",
                    bins,
                    rotationalTemperature_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_rotationalTemperature_3D_pos.xyz",
                    vectorBins,
                    rotationalTemperature_
                );                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_vibrationalTemperature.xy",
                    bins,
                    vibrationalTemperature_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_vibrationalTemperature_3D_pos.xyz",
                    vectorBins,
                    vibrationalTemperature_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_overallTemperature.xy",
                    bins,
                    overallTemperature_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_overallTemperature_3D_pos.xyz",
                    vectorBins,
                    overallTemperature_
                ); 
            }

            // output pressure
            if(outputField_[3])
            {
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_pressureTensor.xyz",
                    bins,
                    pField_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_pressureTensor_3D_pos.xyz",
                    vectorBins,
                    pField_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_p.xy",
                    bins,
                    scalarPressure_
                );
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_p_3D_pos.xyz",
                    vectorBins,
                    scalarPressure_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_stressTensor.xyz",
                    bins,
                    tauField_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_stressTensor_3D_pos.xyz",
                    vectorBins,
                    tauField_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_heatFluxVector.xyz",
                    bins,
                    qField_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_heatFluxVector_3D_pos.xyz",
                    vectorBins,
                    qField_
                ); 
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_internalHeatFluxVector.xyz",
                    bins,
                    qInternalField_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_internalHeatFluxVector_3D_pos.xyz",
                    vectorBins,
                    qInternalField_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_translationalHeatFluxVector.xyz",
                    bins,
                    qTranslationalField_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_translationalHeatFluxVector_3D_pos.xyz",
                    vectorBins,
                    qTranslationalField_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_variableHardSphereMeanFreePath.xy",
                    bins,
                    meanFreePath_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_variableHardSphereMeanFreePath_3D_pos.xyz",
                    vectorBins,
                    meanFreePath_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_Ma.xy",
                    bins,
                    Ma_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_Ma_3D_pos.xyz",
                    vectorBins,
                    Ma_
                );
            }
            
            //permeabilityMeasurements
            if(permeabilityMeasurements_)
            {
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binVolume.xy",
                    bins,
                    binVolume_
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binVolume_3D_pos.xyz",
                    vectorBins,
                    binVolume_
                );
                
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binNumberOfParticles.xy",
                    bins,
                    N_*cloud_.nParticle()
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binNumberOfParticles_3D_pos.xyz",
                    vectorBins,
                    N_*cloud_.nParticle()
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binMass.xy",
                    bins,
                    averageMass_
                    //M_/(N_*cloud_.nParticle())
                );
                
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_binMass_3D_pos.xyz",
                    vectorBins,
                    averageMass_
//                     M_/(N_*cloud_.nParticle())
                );
            }
        }
    }
}

void pdBinsMethod::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //
