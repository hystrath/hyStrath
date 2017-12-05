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

Controls translational temperature in bins in a user-defined zone in a mesh.
Allows for a mixture of cell-by-cell and zone treatment, improved statistics
compared to cell-by-cell and increased resolution compared to zone treatment.

\*---------------------------------------------------------------------------*/

#include "temperatureBerendsenBins.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(temperatureBerendsenBins, 0);

addToRunTimeSelectionTable(dsmcStateController, temperatureBerendsenBins, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
temperatureBerendsenBins::temperatureBerendsenBins
(
    Time& t,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcStateController(t, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    tauT_(readScalar(propsDict_.lookup("tauT"))),
    typeIds_(),
    averagingCounter_(0.0),
    resetFieldsAtOutput_(true)
{

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
            FatalErrorIn("temperatureBerendsenBins::temperatureBerendsenBins()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------

    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh_, propsDict_)
    );
    
    const label& nBins = binModel_->nBins();

    mols_.setSize(nBins, 0.0);
    dsmcMols_.setSize(nBins, 0.0);
    mass_.setSize(nBins, 0.0);
    mcc_.setSize(nBins, 0.0);
    mom_.setSize(nBins, vector::zero);
    UCollected_.setSize(nBins, vector::zero);
    UMean_.setSize(nBins, vector::zero);
    chi_.setSize(nBins, 1.0);
    binTemperatures_.setSize(nBins, 0.0);
    translationalTemperature_.setSize(nBins, 0.0);    
    
    bool readTemperatures = false;
    
    if(propsDict_.found("temperatures"))
    {
        readTemperatures = true;
    }
    
    if(!readTemperatures)
    {
        scalar startTemperature = readScalar(propsDict_.lookup("startTemperature"));
        scalar endTemperature = readScalar(propsDict_.lookup("endTemperature"));
        
        scalar meanTemperature = (startTemperature + endTemperature)*0.5;
        
        translationalTemperature_ = meanTemperature;
        binTemperatures_ = meanTemperature;
    
   
        // set target temperature in bins

        Info << "temperatureBerendsenBins: output target temperatures: " << nl << endl;

        const scalar deltaTBin = (endTemperature - startTemperature)/(nBins - 1);
        
        forAll(binTemperatures_, n)
        {
            binTemperatures_[n] = startTemperature + n*deltaTBin;
            
            Info << "\t bin (#) " << n  << ": T = " << binTemperatures_[n] << endl;
        }
    }
    else
    {
        List<scalar> temperatures (propsDict_.lookup("temperatures"));
        
        if(temperatures.size() != nBins)
        {
            FatalErrorIn("temperatureBerendsenBins::temperatureBerendsenBins()")
                << "bins size = " << nBins << " not equal to temperatures list: "
                << temperatures.size() << nl << "in: "
                << t.system()/"controllersDict"
                << exit(FatalError); 
        }
        
        scalar meanTemperature = 0.0;
        
        forAll(binTemperatures_, n)
        {
            binTemperatures_[n] = temperatures[n];
            meanTemperature += temperatures[n];
            Info << "\t bin (#) " << n  << ": T = " << binTemperatures_[n] << endl;
        }
        
        translationalTemperature_ = meanTemperature/scalar(temperatures.size());
    }
    
    Info << nl << endl;
    
    if(propsDict_.found("resetFieldsAtOutput"))
    {
        resetFieldsAtOutput_= Switch(propsDict_.lookup("resetFieldsAtOutput"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

temperatureBerendsenBins::~temperatureBerendsenBins()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureBerendsenBins::initialConfiguration()
{
}

void temperatureBerendsenBins::calculateProperties()
{
//     timeIndex_++;
    averagingCounter_ += 1.0;
      
    const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
    const labelList& cells = mesh_.cellZones()[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const vector& rI = p->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    const dsmcParcel::constantProperties& constProp 
                                    = cloud_.constProps(p->typeId());
                                    
                    scalar nParticle = cloud_.nParticle();
                    
                    const scalar& RWF = cloud_.getRWF_cell(cellI);
                    
                    nParticle *= RWF;
                    
                    const scalar& mass = constProp.mass()*nParticle;

                    mols_[n] += nParticle;
                    dsmcMols_[n] += 1.0;
                    mass_[n] += mass;
                    mcc_[n] += mass*mag(p->U())*mag(p->U());                   
                    UCollected_[n] += p->U();
                    mom_[n] += mass*p->U();
                }
            }
        }
    }
    
    {
        
        scalarField mass = mass_;
        scalarField mols = mols_;
        scalarField dsmcMols = dsmcMols_;
        scalarField mcc = mcc_;
        vectorField mom = mom_;
        vectorField UCollected = UCollected_;
        
        //- parallel communication

        if(Pstream::parRun())
        {
            forAll(mols, n)
            {
                reduce(mols[n], sumOp<scalar>());
                reduce(dsmcMols[n], sumOp<scalar>());
                reduce(mass[n], sumOp<scalar>());
                reduce(mcc[n], sumOp<scalar>());
                reduce(mom[n], sumOp<vector>());                
                reduce(UCollected[n], sumOp<vector>());
            }
        }

        forAll(mols, n)
        {
            if(mols[n] > 0.0)
            {
                UMean_[n] = UCollected[n]/dsmcMols[n];
                
                translationalTemperature_[n] = (1.0/(3.0*physicoChemical::k.value()))
                    *(
                        ((mcc[n]/(mols[n])))
                        - (
                            (mass[n]/(mols[n])
                            )*mag(UMean_[n])*mag(UMean_[n]))
                    );
            }
        }
        
        if(resetFieldsAtOutput_)
        {
            //- reset fields
            averagingCounter_ = 0.0;
            
            mols_ = 0.0;
            dsmcMols = 0.0;
            mass_ = 0.0;
            mcc_ = 0.0;
            mom_ = vector::zero;
            UCollected_ = vector::zero;
        }
        
        const scalar& deltaTDSMC = mesh_.time().deltaTValue(); // time step 

        Info << " Thermostat (temperatureBerendsenBins) output: " << nl <<  endl;

        chi_ = 1.0;
        
        forAll(translationalTemperature_, n)
        {
            if(translationalTemperature_[n] > 0.0)
            {    
                chi_[n] = sqrt(1.0 + (deltaTDSMC/tauT_)*((binTemperatures_[n]/translationalTemperature_[n]) - 1.0) );

                Info<< "target temperature: " << binTemperatures_[n] 
                    << " measured T: " << translationalTemperature_[n]
                    << " chi: " << chi_[n]
                    << endl;
            }
        }

        Info << nl << endl;
    }
}

void temperatureBerendsenBins::controlParcelsBeforeMove()
{
	
}


void temperatureBerendsenBins::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {

    }
}


void temperatureBerendsenBins::controlParcelsBeforeCollisions()
{
    Info << "temperatureBerendsenBins: control" << endl;

    const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();
            
    const labelList& cells = mesh_.cellZones()[regionId_];

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        
        const List<dsmcParcel*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const vector& rI = p->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(typeIds_, p->typeId()) != -1)
                {
                    if(chi_[n] > 0)
                    {
                        p->U() -= UMean_[n];
                        p->U() *= chi_[n];
                        p->U() += UMean_[n];
                    }
                }
            }
        }
    }
}

void temperatureBerendsenBins::controlParcelsAfterCollisions()
{}

void temperatureBerendsenBins::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

}

void temperatureBerendsenBins::setProperties()
{

}


}
// End namespace Foam

// ************************************************************************* //
