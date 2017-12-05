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

#include "dsmcRayleighBrillouinBinsMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcRayleighBrillouinBinsMethod, 0);

addToRunTimeSelectionTable(dsmcField, dsmcRayleighBrillouinBinsMethod, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcRayleighBrillouinBinsMethod::dsmcRayleighBrillouinBinsMethod
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),
    typeIds_(),

    averagingCounter_(0.0),
    mols_(),
    rhoN_(),
    outputField_(4, true),
    averagingAcrossManyRuns_(false)
    
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcRayleighBrillouinBinsMethod::dsmcRayleighBrillouinBinsMethod()")
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
    
    // ---------------------------------------------------

    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );
    
    const label& nBins = binModel_->nBins();

    mols_.setSize(nBins, 0.0);
    
    rhoN_.setSize(nBins, 0.0);
    
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcRayleighBrillouinBinsMethod::~dsmcRayleighBrillouinBinsMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcRayleighBrillouinBinsMethod::readIn()
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

    dict.readIfPresent("mols", mols_);
    dict.readIfPresent("averagingCounter", averagingCounter_);
}

void dsmcRayleighBrillouinBinsMethod::writeOut()
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
        
        dict.add("mols", mols_);
        dict.add("averagingCounter", averagingCounter_);
        
        IOstream::streamFormat fmt = time_.time().writeFormat();

        IOstream::versionNumber ver = time_.time().writeVersion();

        IOstream::compressionType cmp = time_.time().writeCompression();
    
        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

void dsmcRayleighBrillouinBinsMethod::createField()
{  
}


void dsmcRayleighBrillouinBinsMethod::calculateField()
{
    if(time_.averagingTime())
    {
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
                
                label iD = findIndex(typeIds_, p->typeId());

                const vector& rI = p->position();

                label n = binModel_->isPointWithinBin(rI, cellI);

                if(n != -1)
                {
                    if(iD != -1)
                    {
                        if(cloud_.axisymmetric())
                        {
                            const point& cC = cloud_.mesh().cellCentres()[cellI];
                            scalar radius = cC.y();
                            
                            scalar RWF = 1.0 + cloud_.maxRWF()*(radius/cloud_.radialExtent());
                            
                            mols_[n] += cloud_.nParticle()*RWF;
                        }
                        else
                        {
                            mols_[n] += cloud_.nParticle();
                        }
                    }
                }
            }
        }

        scalarField mols = mols_;
        
        //- parallel communication

        if(Pstream::parRun())
        {
            forAll(mols, n)
            {
                reduce(mols[n], sumOp<scalar>());
            }
        }

        forAll(mols, n)
        {
            scalar volume = binModel_->binVolume(n);

            rhoN_[n] = (mols[n])/(averagingCounter_*volume);
        }

        if(time_.resetFieldsAtOutput())
        {
            //- reset fields
            averagingCounter_ = 0.0;
            
            mols_ = 0.0;    
        }
        
        if(averagingAcrossManyRuns_)
        {
            writeOut();
        }
    }
}

void dsmcRayleighBrillouinBinsMethod::writeField()
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
            }
        }
    }
}

void dsmcRayleighBrillouinBinsMethod::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}


} // End namespace Foam

// ************************************************************************* //
