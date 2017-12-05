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

Measures translational, rotational, vibrational and overall temperature in a user-defined zone.

Results are written to timeDirectory/uniform as text files.

Time averaging is switched on or off with the resetAtOutput option in fieldPropertiesDict.

\*---------------------------------------------------------------------------*/

#include "dsmcVibrationalDegreesOfFreedomZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVibrationalDegreesOfFreedomZone, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVibrationalDegreesOfFreedomZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVibrationalDegreesOfFreedomZone::dsmcVibrationalDegreesOfFreedomZone
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
	vibrationalETotal_(),
    nParcels_(),
    vibT_(),
    vDof_(),
    vibTxvDof_(0.0),
    vDoF_(0.0),
    vDoFField_(time_.totalNAvSteps()+1, 0.0)
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
            FatalErrorIn("dsmcVibrationalDegreesOfFreedomZone::dsmcVibrationalDegreesOfFreedomZone()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << t.system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }
    // ---------------------------------------------------


    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("dsmcVibrationalDegreesOfFreedomZone::dsmcVibrationalDegreesOfFreedomZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
    vibrationalETotal_.setSize(typeIds_.size());
    
    nParcels_.setSize(typeIds_.size());
    
    vibT_.setSize(typeIds_.size());
    
    vDof_.setSize(typeIds_.size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVibrationalDegreesOfFreedomZone::~dsmcVibrationalDegreesOfFreedomZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- initial conditions
void dsmcVibrationalDegreesOfFreedomZone::createField()
{
    Info << "Initialising dsmcVibrationalDegreesOfFreedomZone field" << endl;
    
    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize(cloud_.constProps(typeIds_[i]).vibrationalDegreesOfFreedom(), 0.0);
    }
}


void dsmcVibrationalDegreesOfFreedomZone::calculateField()
{
    if(time_.samplingTime())
    {
        const List< DynamicList<dsmcParcel*> >& cellOccupancy
            = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        forAll(cells, c)
        {
            const label& cellI = cells[c];
            const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cellI];
            
            forAll(vibrationalETotal_, iD)
            {
                forAll(parcelsInCell, pIC)
                {
                    dsmcParcel* p = parcelsInCell[pIC];
    
                    if(p->typeId() == typeIds_[iD])
                    {
//                         vibrationalETotal_[iD] += p->vibLevel()*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV();
                        forAll(vibrationalETotal_[iD], m)
                        {
                             vibrationalETotal_[iD][m] += p->vibLevel()[m]*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV()[m];
                        }
                        nParcels_[iD] += 1.0;
                    }
                } 
            }  
        }   
    }

    if(time_.averagingTime())
    {

        const scalar& timeIndex = time_.averagingTimeIndex();
        
        List<scalarList> degreesOfFreedomMode;
        List<scalarList> vibTMode;
        
        
        degreesOfFreedomMode.setSize(typeIds_.size());
        vibTMode.setSize(typeIds_.size());
                
        forAll(vibrationalETotal_, iD)
        {
//             if(vibrationalETotal_[iD] > VSMALL && nParcels_[iD] > VSMALL)
//             {
//                 const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
//                 
//                 scalar vibrationalEMean = (vibrationalETotal_[iD]/nParcels_[iD]);
//                 scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//                 
//                 vibT_[iD] = thetaV / log(1.0 + (1.0/iMean));
//                 vDof_[iD] = (2.0*thetaV/vibT_[iD]) / (exp(thetaV/vibT_[iD]) - 1.0);
//                 
//                 vDoF_ += vDof_[iD];
//             } 
            
            forAll(vibrationalETotal_[iD], v)
            {
                if(vibrationalETotal_[iD][v] > VSMALL && nParcels_[iD] > VSMALL)
                {        
                    scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV()[v];
                    
                    scalar vibrationalEMean = (vibrationalETotal_[iD][v]/nParcels_[iD]);
                    
                    scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
                    
                    vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                    degreesOfFreedomMode[iD][v] = (2.0*thetaV/vibTMode[iD][v]) / (exp(thetaV/vibTMode[iD][v]) - 1.0);
                }
            }
            
            forAll(degreesOfFreedomMode[iD], v)
            {
                vDoF_ += degreesOfFreedomMode[iD][v];
            }
        }
        
        vDoFField_[timeIndex] = vDoF_;
        
        vDoF_ = 0.0; //reset this so it doesn't keep accumulating, which would bias the overallT towards vibrationalT
        
        if(Pstream::parRun())
        {
            reduce(vDoFField_[timeIndex], sumOp<scalar>());
        }
    }
}

//- write field
void dsmcVibrationalDegreesOfFreedomZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");
    
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        const scalarField& timeField = time_.averagingTimes();    

        writeTimeData(timePath, fieldName_+"_"+regionName_+"_vibrationalDegreesOfFreedom", timeField, vDoFField_);

        
        //- reset
        if(time_.resetFieldsAtOutput())
        {
            nParcels_ = scalar(0.0);
            vibT_ = scalar(0.0);
            vDof_ = scalar(0.0);
            vibTxvDof_ = scalar(0.0);
            
            forAll(vibrationalETotal_, m)
            {
                vibrationalETotal_[m] = scalar(0.0);
            }
        }
    }
}


void dsmcVibrationalDegreesOfFreedomZone::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

}

} // End namespace Foam

// ************************************************************************* //
