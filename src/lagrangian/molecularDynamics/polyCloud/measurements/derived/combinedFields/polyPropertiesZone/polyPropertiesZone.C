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

#include "polyPropertiesZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyPropertiesZone, 0);
addToRunTimeSelectionTable(polyField, polyPropertiesZone, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyPropertiesZone::polyPropertiesZone
(
    Time& t,
    const polyMesh& mesh,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyField(t, mesh, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fields_(t, mesh, "dummy"),
    fieldName_(propsDict_.lookup("fieldName")),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    totalVolume_(0.0),
    mols_(0.0),
    mass_(0.0),
    mom_(vector::zero),
    kE_(0.0),
    velocityB_(vector::zero),
    virialTensor_(tensor::zero),
    kineticTensor_(tensor::zero),
    kineticTensor2_(tensor::zero),
    dof_(0.0),
    keSum_(0.0),
    peSum_(0.0),
    angularKeSum_(0.0),
    molIds_(),
    timeIndex_(0),
	molField_(1, 0.0),
	densityField_(1, 0.0),
	massDensityField_(1, 0.0),
	velocityFieldA_(1, vector::zero),
	velocityFieldB_(1, vector::zero),
	momentumField_(1, vector::zero),
	temperatureField_(1, 0.0),
	stressField_(1, tensor::zero),
	pressureField_(1, 0.0),
	stressField2_(1, tensor::zero),
	pressureField2_(1, 0.0),
	kEField_(1, 0.0),
	pEField_(1, 0.0),
	angularkEField_(1, 0.0),
	tEField_(1, 0.0),
    outputField_(5, true),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyPropertiesZone::polyPropertiesZone()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
    bool readFromStore = true;
    
    if (propsDict_.found("readFromStorage"))
    {
        readFromStore = Switch(propsDict_.lookup("readFromStorage"));
    }    

    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));  

	if (!resetAtOutput_ && readFromStore)
    {
        Info << " Averaging across many runs. Reading from dictionary:" << endl;

        pathName_ = time_.time().path()/"storage";
        nameFile_ = "propertiesZoneData_"+fieldName_;

        if( !isDir(pathName_) )
        {
            mkDir(pathName_);

            Info << nl << "Storage not found!"  << nl << endl;
            Info << ".... creating"  << nl << endl;
        }

        bool fileFound = readFromStorage();

        if(!fileFound)
        {
            Info << nl << "File not found: " << nameFile_ << nl << endl;
            Info << ".... creating"  << nl << endl;
            writeToStorage();

            Info << "setting properties to default values. " << endl;
        }
        else
        {
            Info << "Reading from storage, e.g. noAvTimeSteps = " << nAvTimeSteps_ << endl;            
//             Pout<< "Properties read-in are: mols = " << mols_ << ", mass = " << mass_
//                 << ", averagingTime = " << nAvTimeSteps_
//                 << endl;
        }
       
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

    // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();


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

        //propertyNames.shrink();

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

        if(findIndex(propertyNames, "energy") == -1)
        {
            outputField_[4] = false;
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
                (propertyName != "pressure") &&
                (propertyName != "energy")
            )
            {    
                FatalErrorIn("atomisticBinsMethod::atomisticBinsMethod()")
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << time_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);            
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyPropertiesZone::~polyPropertiesZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyPropertiesZone::createField()
{}

void polyPropertiesZone::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    scalar mols = 0.0;
    scalar mass = 0.0;
    vector mom = vector::zero;
    vector angularSpeed = vector::zero;

    forAll(mesh_.cellZones()[regionId_], c) 
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];
    
        forAll(molsInCell, m)
        {
            polyMolecule* molI = molsInCell[m];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                mols += 1.0;
                const scalar& massI = molCloud_.cP().mass(molI->id());
                mass += massI;
                mom += massI*molI->v();

                const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                // angular speed 
                const vector& molOmega(inv(molMoI) & molI->pi());

                angularSpeed += molOmega;
            }
        }
    }


    // - parallel processing
    if(Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());            
        reduce(mass, sumOp<scalar>());
        reduce(mom, sumOp<vector>());
        reduce(angularSpeed, sumOp<vector>());
    }
    
    vector velocity = vector::zero;
    
    if(mass > 0.0)
    {
        velocity = mom/mass;
    }
    
    velocityB_ += velocity;    
    
    vector angularVelocity = vector::zero;
    
    if(mols > 0.0)
    {
        angularVelocity = angularSpeed/mols;
    }

    mols = 0.0;
    mass = 0.0;
    mom = vector::zero;
    scalar kE = 0.0;
    tensor virialTensor = tensor::zero;
    tensor kineticTensor = tensor::zero;
    tensor kineticTensor2 = tensor::zero;
    scalar dof = 0.0;
    scalar keSum = 0.0;
    scalar peSum = 0.0;
    scalar angularKeSum = 0.0;

    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            if(findIndex(molIds_, molI->id()) != -1)
            {
                mols += 1.0;
                const scalar& massI = molCloud_.cP().mass(molI->id());

                mass += massI;
                mom += massI*molI->v();

                dof += molCloud_.cP().degreesOfFreedom(molI->id());

                kE += 0.5*massI*magSqr(molI->v() - velocity);

                keSum += 0.5*massI*magSqr(molI->v());
                peSum += molI->potentialEnergy();

                const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                // angular speed 
                const vector& molOmega(inv(molMoI) & molI->pi());
                angularKeSum += 0.5*(molOmega & molMoI & molOmega);


                kineticTensor += ( massI*(molI->v() - velocity)*(molI->v() - velocity) ) 
                                    + ( ((molOmega - angularVelocity) & molMoI) * (molOmega-angularVelocity) );

                kineticTensor2 += ( massI*(molI->v() - velocity)*(molI->v() - velocity) );

                virialTensor += 0.5*molI->rf();
            }
        }
    }

    mols_ += mols;
    mass_ += mass;
    mom_ += mom;
    kE_ += kE;
    virialTensor_ += virialTensor;
    kineticTensor_ += kineticTensor;
    kineticTensor2_ += kineticTensor2;
    dof_ += dof;
    keSum_ += keSum;
    peSum_ += peSum;
    angularKeSum_ += angularKeSum;

    if(time_.outputTime()) 
    {
        //- parallel communication
        scalar mols = mols_;
        scalar mass = mass_;
        vector mom = mom_;
        vector velocityB = velocityB_;
        scalar kE = kE_;
        scalar dof = dof_;
        tensor virialTensor = virialTensor_;
        tensor kineticTensor = kineticTensor_;
        tensor kineticTensor2 = kineticTensor2_;
        scalar keSum = keSum_;
        scalar peSum = peSum_;
        scalar angularKeSum = angularKeSum_;

        if(Pstream::parRun())
        {
            reduce(mols, sumOp<scalar>());
            reduce(mass, sumOp<scalar>());
            reduce(mom, sumOp<vector>());
            reduce(kE, sumOp<scalar>());
            reduce(dof, sumOp<scalar>());
            reduce(kineticTensor, sumOp<tensor>());
            reduce(kineticTensor2, sumOp<tensor>());
            reduce(virialTensor, sumOp<tensor>());
            reduce(keSum, sumOp<scalar>());
            reduce(peSum, sumOp<scalar>());
            reduce(angularKeSum, sumOp<scalar>());
        }

        const scalar& nAvTimeSteps = nAvTimeSteps_; 

        // density
        
        molField_[timeIndex_] = mols/nAvTimeSteps;
        densityField_[timeIndex_] = mols/(totalVolume_*nAvTimeSteps);
        massDensityField_[timeIndex_] = mass/(totalVolume_*nAvTimeSteps);


        // velocity 
        vector velocityA = vector::zero;

        if(mass > 0.0)
        {
            velocityA = mom/mass;
        }

        velocityFieldA_[timeIndex_] = velocityA;
        velocityFieldB_[timeIndex_] = velocityB/nAvTimeSteps;

        vector momentum = vector::zero;

        if(mols > 0.0)
        {
            momentum = mom/mols;
        }

        momentumField_[timeIndex_] = momentum;


        // temperature
        const scalar& kB = molCloud_.redUnits().kB();

        scalar zoneTemperature = 0.0;

        if(dof > 0.0)
        {
            zoneTemperature = (2.0*(kE+angularKeSum))/(kB*dof);
        }

        temperatureField_[timeIndex_] = zoneTemperature;

        // pressure

        scalar zonePressure = 0.0;
        tensor zoneStress = tensor::zero;

        scalar zonePressure2 = 0.0;
        tensor zoneStress2 = tensor::zero;

        if(dof > 0.0)
        {
            zonePressure = tr( (3.0*mols*kineticTensor/dof) + virialTensor)
                                    /( 3.0*totalVolume_*nAvTimeSteps );

            zoneStress = ((3.0*mols*kineticTensor/dof) + virialTensor)
                                    /(totalVolume_*nAvTimeSteps);

            zonePressure2 = tr( (3.0*mols*kineticTensor2/dof) + virialTensor)
                                    /( 3.0*totalVolume_*nAvTimeSteps );

            zoneStress2 = ((3.0*mols*kineticTensor2/dof) + virialTensor)
                                    /(totalVolume_*nAvTimeSteps);
        }

        stressField_[timeIndex_] = zoneStress;

        pressureField_[timeIndex_] = zonePressure;

        stressField2_[timeIndex_] = zoneStress2;

        pressureField2_[timeIndex_] = zonePressure2;

        // energies

        scalar kineticEnergy = 0.0;
        scalar potentialEnergy = 0.0;
        scalar angularKineticEnergy = 0.0;

        if(mols > 0.0)
        {
            kineticEnergy = keSum/mols;
            potentialEnergy = peSum/mols;
            angularKineticEnergy = angularKeSum_/mols;
        }

        kEField_[timeIndex_] = kineticEnergy;
        pEField_[timeIndex_] = potentialEnergy;
        angularkEField_[timeIndex_] = angularKineticEnergy;
        tEField_[timeIndex_] = kineticEnergy + potentialEnergy + angularKineticEnergy;


        if(resetAtOutput_)
        {
            nAvTimeSteps_ = 0.0;
            mols_ = 0.0;
            mass_ = 0.0;
            mom_ = vector::zero;
            velocityB_ = vector::zero;
            kE_ = 0.0;
            dof_ = 0.0;
            virialTensor_ = tensor::zero;
            kineticTensor_ = tensor::zero;
            kineticTensor2_ = tensor::zero;
            keSum_ = 0.0;
            peSum_ = 0.0;
            angularKeSum_ = 0.0;
        }
        else 
        {
            writeToStorage();
        }
    }
}

void polyPropertiesZone::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
        file << nAvTimeSteps_ << endl;
        file << mols_ << endl;
        file << mass_ << endl;
        file << mom_ << endl; 
        file << kE_ << endl;        
        file << velocityB_ << endl;
        file << virialTensor_ << endl;
        file << kineticTensor_ << endl;
        file << kineticTensor2_ << endl;
        file << dof_ << endl;        
        file << keSum_ << endl;
        file << peSum_ << endl;
        file << angularKeSum_ << endl;

              
    }
    else
    {
        FatalErrorIn("void polyPropertiesZone::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool polyPropertiesZone::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar nAvTimeSteps;
        scalar mols;
        scalar mass;
        vector mom;
        scalar kE;        
        vector velocityB;
        tensor virialTensor;
        tensor kineticTensor;
        tensor kineticTensor2;
        scalar dof;        
        scalar keSum;
        scalar peSum;
        scalar angularKeSum;

        file >> nAvTimeSteps;
        file >> mols;
        file >> mass;
        file >> mom;
        file >> kE;
        file >> velocityB;
        file >> virialTensor;
        file >> kineticTensor;
        file >> kineticTensor2;
        file >> dof;        
        file >> keSum;
        file >> peSum;
        file >> angularKeSum;

        nAvTimeSteps_ = nAvTimeSteps;
        mols_ = mols;
        mass_ = mass;
        mom_ = mom;
        kE_ = kE;
        velocityB_ = velocityB;
        virialTensor_ = virialTensor;        
        kineticTensor_ = kineticTensor;
        kineticTensor2_ = kineticTensor2;
        dof_ = dof;        
        keSum_ = keSum;
        peSum_ = peSum;
        angularKeSum_ = angularKeSum;
    }

    return goodFile;    
}


void polyPropertiesZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());

            // output densities
            if(outputField_[0])
            {

                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_N.xy",
                    timeField,
                    molField_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_rhoN.xy",
                    timeField,
                    densityField_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_rhoM.xy",
                    timeField,
                    massDensityField_,
                    true
                );
            }

            // output velocities
            if(outputField_[1])
            {

            	writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_U_CAM.xyz",
                    timeField,
                    velocityFieldA_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_U_SAM.xyz",
                    timeField,
                    velocityFieldB_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_mom.xyz",
                    timeField,
                    momentumField_,
                    true
                );
            }
            // output temperature
            if(outputField_[2])
            {

            	writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_T.xy",
                    timeField,
                    temperatureField_,
                    true
                );
            }

            // output pressure
            if(outputField_[3])
            {

            	writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_stress.xyz",
                    timeField,
                    stressField_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_p.xy",
                    timeField,
                    pressureField_,
                    true
                );

                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_stress2.xyz",
                    timeField,
                    stressField2_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_p2.xy",
                    timeField,
                    pressureField2_,
                    true
                );
            }

            // output energies
            if(outputField_[4])
            {

            	writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_kE.xy",
                    timeField,
                    kEField_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_pE.xy",
                    timeField,
                    pEField_,
                    true
                );
    
                writeTimeData
                (
                    casePath_,
                    "zone_"+regionName_+"_"+fieldName_+"_tE.xy",
                    timeField,
                    tEField_,
                    true
                );
            }

            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                // output densities
                if(outputField_[0])
                {

                	writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_N_SI.xy",
                        timeField*rU.refTime(),
                        molField_,
                        true
                    );
    
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_rhoN_SI.xy",
                        timeField*rU.refTime(),
                        densityField_*rU.refNumberDensity(),
                        true
                    );
    
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_rhoM_SI.xy",
                        timeField*rU.refTime(),
                        massDensityField_*rU.refMassDensity(),
                        true
                    );
                }

                // output velocities
                if(outputField_[1])
                {

                	writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_U_CAM_SI.xyz",
                        timeField*rU.refTime(),
                        velocityFieldA_*rU.refVelocity(),
                        true
                    );
        
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_U_SAM_SI.xyz",
                        timeField*rU.refTime(),
                        velocityFieldB_*rU.refVelocity(),
                        true
                    );
                }
                // output temperature
                if(outputField_[2])
                {

                	writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_T_SI.xy",
                        timeField*rU.refTime(),
                        temperatureField_*rU.refTemp(),
                        true
                    );
                }

                // output pressure
                if(outputField_[3])
                {
                	writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_stress_SI.xyz",
                        timeField*rU.refTime(),
                        stressField_*rU.refPressure(),
                        true
                    );
        
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_p_SI.xy",
                        timeField*rU.refTime(),
                        pressureField_*rU.refPressure(),
                        true
                    );
                }
                // output energies
                if(outputField_[4])
                {

                	writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_kE_SI.xy",
                        timeField*rU.refTime(),
                        kEField_*rU.refEnergy(),
                        true
                    );
        
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_pE_SI.xy",
                        timeField*rU.refTime(),
                        pEField_*rU.refEnergy(),
                        true
                    );
        
                    writeTimeData
                    (
                        casePath_,
                        "zone_"+regionName_+"_"+fieldName_+"_tE_SI.xy",
                        timeField*rU.refTime(),
                        tEField_*rU.refEnergy(),
                        true
                    );
                }
            }
        }
    }
}

void polyPropertiesZone::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyPropertiesZone::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyPropertiesZone::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
