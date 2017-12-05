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

#include "polyBinsMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyBinsMethod, 0);
addToRunTimeSelectionTable(polyField, polyBinsMethod, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyBinsMethod::polyBinsMethod
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
    binModel_(),
    regionName_(propsDict_.lookup("zoneName")),
    regionId_(-1),
    fieldName_(propsDict_.lookup("fieldName")),

    molIds_(),

    mols_(),
    mass_(),
    mom_(),
    velocityB_(),
    kE_(),
    angularKeSum_(),
    dof_(),
    kineticTensor_(),
    virialTensor_(),

    molsV_(),
    massV_(),
    momV_(),
    velocity_(),
    angularSpeed_(),
    angularVelocity_(),

    N_(),
    rhoN_(),
    rhoM_(),
    USAM_(),
    UCAM_(),
    T_(),
    stress_(),
    p_(),
    outputField_(4, true),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
 
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("polyBinsMethod::polyBinsMethod()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
    // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    
    // create bin model
    binModel_ = autoPtr<binModel>
    (
        binModel::New(mesh, propsDict_)
    );

    const label& nBins = binModel_->nBins();

    mols_.setSize(nBins, 0.0);
    mass_.setSize(nBins, 0.0);
    mom_.setSize(nBins, vector::zero);
    velocityB_.setSize(nBins, vector::zero);
    kE_.setSize(nBins, 0.0);
    angularKeSum_.setSize(nBins, 0.0);
    dof_.setSize(nBins, 0.0);
    kineticTensor_.setSize(nBins, tensor::zero);
    virialTensor_.setSize(nBins, tensor::zero);

    molsV_.setSize(nBins, 0.0);
    massV_.setSize(nBins, 0.0);
    momV_.setSize(nBins, vector::zero);
    velocity_.setSize(nBins, vector::zero);
    angularSpeed_.setSize(nBins, vector::zero);
    angularVelocity_.setSize(nBins, vector::zero);

    N_.setSize(nBins, 0.0);
    rhoN_.setSize(nBins, 0.0);
    rhoM_.setSize(nBins, 0.0);
    USAM_.setSize(nBins, vector::zero);
    UCAM_.setSize(nBins, vector::zero);
    T_.setSize(nBins, 0.0);
    stress_.setSize(nBins, tensor::zero);
    p_.setSize(nBins, 0.0);
    pVir_.setSize(nBins, 0.0);
    pKin_.setSize(nBins, 0.0);


  
    {
        Info << nl << "Storage..." << endl;
        
        bool resetStorage = false;
        
        if (propsDict_.found("resetStorage"))
        {
            resetStorage = Switch(propsDict_.lookup("resetStorage"));
        }    
        
        if(resetStorage)
        {
            Info<< "WARNING: storage will be reset."
                << " This is not good if you would like to average over many runs. "
                << endl;
        }
        else
        {
            Info << "WARNING: storage will NOT be reset."
                << " This is good if you would like to average over many runs. "
                << " This is NOT good if you have been testing your simulation a number of times "
                << " Delete your storage directory before moving to important runs"
                << " or type resetStorage = yes, for just the first simulation run."
                << endl;            
        }
        
        resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));    
        
        
        // stored data activation in dictionary        

        pathName_ = time_.time().path()/"storage";
        nameFile_ = "binsData_"+fieldName_;

        if( !isDir(pathName_) )
        {
            mkDir(pathName_);

            Info << nl << "Storage not found!"  << nl << endl;
            Info << ".... creating"  << nl << endl;
        }

        if(!resetStorage)
        {
            readFromStorage();
        }

        IFstream file(pathName_/nameFile_);

        bool foundFile = file.good();
        
        if(!foundFile)
        {
            Info << nl << "File not found: " << nameFile_ << nl << endl;
            Info << ".... creating"  << nl << endl;
            writeToStorage();

            Info << "setting properties to default values. " << endl;
        }
        else
        {
            Info << "Reading from storage, e.g. noAvTimeSteps = " << nAvTimeSteps_ << endl;
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
                FatalErrorIn("polyBinsMethod::polyBinsMethod()")
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << time_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);            
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyBinsMethod::~polyBinsMethod()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyBinsMethod::createField()
{
    Info << "Initialising polyBinsMethod fields" << endl;

    const List< DynamicList<polyMolecule*> >& cellOccupancy
        = molCloud_.cellOccupancy();

    const labelList& cells = mesh_.cellZones()[regionId_];

    scalarField mols(mass_.size(), 0.0);
    scalarField mass(mass_.size(), 0.0);
    vectorField mom(mom_.size(), vector::zero);
    vectorField angularSpeed(mom_.size(), vector::zero);

    forAll(cells, c)
    {
        const label& cellI = cells[c];
        const List<polyMolecule*>& molsInCell = cellOccupancy[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
//                     const polyMolecule::constantProperties& constProp 
//                                 = molCloud_.constProps(molI->id());

                    mols[n] += 1.0;

                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    
                    mass[n] += massI;

                    mom[n] += massI*molI->v();

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularSpeed[n] += molOmega;
                }
            }
        }
    }


    // parallel processing
    if(Pstream::parRun())
    {
        // sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << mols << mass << mom << angularSpeed;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalarField molsProc;
                scalarField massProc;
                vectorField momProc;
                vectorField angularSpeedProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> molsProc >> massProc >> momProc >> angularSpeedProc;
                }

                mols += molsProc;
                mass += massProc;
                mom += momProc;
                angularSpeed += angularSpeedProc;
            }
        }
    }

    forAll(velocity_, n)
    {
        if(mass[n] > 0.0)
        {
            angularVelocity_[n] = angularSpeed[n]/mols[n];
            velocity_[n] = mom[n]/mass[n];
        }
    }
}


void polyBinsMethod::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    molsV_[n] += 1.0;
                    
                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    
                    massV_[n] += massI;

                    momV_[n] += massI*molI->v();

                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularSpeed_[n] += molOmega;
                }
            }
        }
    }

    if(time_.outputTime()) 
    {
        scalarField mols = molsV_;
        scalarField mass = massV_;
        vectorField mom = momV_;
        vectorField angularSpeed = angularSpeed_;

        //- parallel communication
        if(Pstream::parRun())
        {
            // sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << mols << mass << mom << angularSpeed;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField molsProc;
                    scalarField massProc;
                    vectorField momProc;
                    vectorField angularSpeedProc;
    
                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> molsProc >> massProc >> momProc >> angularSpeedProc;
                    }
    
                    mols += molsProc;
                    mass += massProc;
                    mom += momProc;
                    angularSpeed += angularSpeedProc;
                }
            }
        }

        velocity_ = vector::zero;
        angularVelocity_ = vector::zero;

        forAll(velocity_, n)
        {
            if(mass[n] > 0.0)
            {
                velocity_[n] = mom[n]/mass[n];
                angularVelocity_[n] = angularSpeed[n]/mols[n];
            }
        }
    }


    scalarField mass(mass_.size(), 0.0);
    vectorField mom(mom_.size(), vector::zero);

    forAll(mesh_.cellZones()[regionId_], c) 
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            label n = binModel_->isPointWithinBin(rI, cellI);

            if(n != -1)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
//                     const polyMolecule::constantProperties& constProp 
//                                             = molCloud_.constProps(molI->id());

                    const scalar& massI = molCloud_.cP().mass(molI->id());
                    mols_[n] += 1.0;
                    mass[n] += massI;
                    mom[n] += massI*molI->v();

                    kE_[n] += 0.5*massI*magSqr(molI->v() - velocity_[n]);


                    dof_[n] += molCloud_.cP().degreesOfFreedom(molI->id());


                    const diagTensor& molMoI(molCloud_.cP().momentOfInertia(molI->id()));

                    // angular speed 
                    const vector& molOmega(inv(molMoI) & molI->pi());

                    angularKeSum_[n] += 0.5*(molOmega & molMoI & molOmega);

                    kineticTensor_[n] += ( massI*(molI->v() - velocity_[n])*(molI->v() - velocity_[n]) ) 
                                            + ( ((molOmega - angularVelocity_[n]) & molMoI)
                                            *(molOmega-angularVelocity_[n]) 
                                            );

                    virialTensor_[n] += 0.5*molI->rf();
                }
            }
        }
    }

    mass_ += mass;
    mom_ += mom;

    if(Pstream::parRun())
    {
        //-sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << mass << mom;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalarField massProc;
                vectorField momProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour  >> massProc >> momProc;
                }
    
                mass += massProc;
                mom += momProc;
            }
        }
    }

    forAll(velocity_, n)
    {
        if(mass[n] > 0.0)
        {
            velocityB_[n] += mom[n]/mass[n];
        }
    }
 
    if(time_.outputTime())
    {
        scalarField mass = mass_;
        scalarField mols = mols_;
        vectorField mom = mom_;
        scalarField kE = kE_;
        scalarField angularKeSum = angularKeSum_;
        scalarField dof = dof_;
        tensorField kineticTensor = kineticTensor_;
        tensorField virialTensor = virialTensor_;

        //- parallel communication
        if(Pstream::parRun())
        {
            //-sending
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);

                        toNeighbour << mols << mass << mom 
                                    << kE << angularKeSum << dof << kineticTensor 
                                    << virialTensor;
                    }
                }
            }
        
            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if(p != Pstream::myProcNo())
                {
                    scalarField molsProc;
                    scalarField massProc;
                    vectorField momProc;
                    scalarField kEProc;
                    scalarField angularKeSumProc;
                    scalarField dofProc;
                    tensorField kineticTensorProc;
                    tensorField virialTensorProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour  >> molsProc >> massProc >> momProc 
                                    >> kEProc >> angularKeSumProc >> dofProc >> kineticTensorProc
                                    >> virialTensorProc;
                    }
        
                    mols += molsProc;
                    mass += massProc;
                    mom += momProc;
                    kE += kEProc;
                    angularKeSum += angularKeSumProc;
                    dof += dofProc;
                    kineticTensor += kineticTensorProc;
                    virialTensor += virialTensorProc;
                }
            }
        }

        const scalar& kB = molCloud_.redUnits().kB();

        const scalar& nAvTimeSteps = nAvTimeSteps_;
        
        forAll(mols, n)
        {
            scalar volume = binModel_->binVolume(n);

            N_[n] = mols[n]/nAvTimeSteps;
            rhoN_[n] = mols[n]/(nAvTimeSteps*volume);
            rhoM_[n] = mass[n]/(nAvTimeSteps*volume);

            USAM_[n] = velocityB_[n]/nAvTimeSteps;

            UCAM_[n] = vector::zero;

            if(mass[n] > 0.0)
            {
                UCAM_[n] = mom[n]/mass[n];
            }

            if(dof[n] > 0.0)
            {
                T_[n] = (2.0*(kE[n]+angularKeSum[n]))/(kB*dof[n]);
            }

            if(dof[n] > 0.0)
            {
                p_[n] = tr( (3.0*mols[n]*kineticTensor[n]/dof[n]) + virialTensor[n])
                                        /( 3.0*volume*nAvTimeSteps );

                pVir_[n] = tr(virialTensor[n])
                                        /( 3.0*volume*nAvTimeSteps );
                pKin_[n] = tr( (3.0*mols[n]*kineticTensor[n]/dof[n]) )
                                        /( 3.0*volume*nAvTimeSteps );

                stress_[n] = ((3.0*mols[n]*kineticTensor[n]/dof[n]) + virialTensor[n])
                                        /( volume*nAvTimeSteps );
            }
        }

        if(resetAtOutput_)
        {
            //- reset fields
            nAvTimeSteps_ = 0.0;
            mols_ = 0.0;
            mass_ = 0.0;
            mom_ = vector::zero;
            velocityB_ = vector::zero;
            kE_ = 0.0;
            angularKeSum_ = 0.0;
            dof_ = 0.0;
            kineticTensor_ = tensor::zero;
            virialTensor_ = tensor::zero;
        }
        else //if(averagingAcrossManyRuns_)
        {
            writeToStorage();
        }
    }
}

void polyBinsMethod::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
        file << nAvTimeSteps_ << endl;
        file << mols_ << endl;
        file << mass_ << endl;
        file << mom_ << endl; 
        file << velocityB_ << endl;
        file << kE_ << endl;
        file << angularKeSum_ << endl;
        file << dof_ << endl;
        file << kineticTensor_ << endl;
        file << virialTensor_ << endl;
    }
    else
    {
        FatalErrorIn("void polyBinsMethod::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool polyBinsMethod::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar nAvTimeSteps;
        scalarField mols;
        scalarField mass;
        vectorField mom;
        vectorField velocityB;
        scalarField kE;
        scalarField angularKeSum;
        scalarField dof;
        tensorField kineticTensor;
        tensorField virialTensor;

        file >> nAvTimeSteps;
        file >> mols;
        file >> mass;
        file >> mom;
        file >> velocityB;
        file >> kE;
        file >> angularKeSum;
        file >> dof;
        file >> kineticTensor;
        file >> virialTensor;

        nAvTimeSteps_ = nAvTimeSteps;
        mols_ = mols;
        mass_ = mass;
        mom_ = mom;
        velocityB_ = velocityB;
        kE_ = kE;
        angularKeSum_ = angularKeSum;
        dof_ = dof;
        kineticTensor_ = kineticTensor;
        virialTensor_ = virialTensor;
    }

    return goodFile;
}

void polyBinsMethod::writeField()
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
                    USAM_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_SAM_3D_pos.xyz",
                    vectorBins,
                    USAM_
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
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_T.xy",
                    bins,
                    T_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_T_3D_pos.xy",
                    vectorBins,
                    T_
                );  
            }

            // output pressure
            if(outputField_[3])
            {
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_stress.xyz",
                    bins,
                    stress_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_stress_3D_pos.xyz",
                    vectorBins,
                    stress_
                );  
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_p.xy",
                    bins,
                    p_
                );
    
                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_pVir.xy",
                    bins,
                    pVir_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_pKin.xy",
                    bins,
                    pKin_
                );

                writeTimeData
                (
                    timePath_,
                    "bins_OneDim_"+regionName_+"_"+fieldName_+"_p_3D_pos.xy",
                    vectorBins,
                    p_
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
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_N_SI.xy",
                        bins*rU.refLength(),
                        N_
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_N_3D_pos_SI.xy",
                        vectorBins*rU.refLength(),
                        N_
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoN_SI.xy",
                        bins*rU.refLength(),
                        rhoN_*rU.refNumberDensity()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoN_3D_pos_SI.xy",
                        vectorBins*rU.refLength(),
                        rhoN_*rU.refNumberDensity()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoM_SI.xy",
                        bins*rU.refLength(),
                        rhoM_*rU.refMassDensity()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_rhoM_3D_pos_SI.xy",
                        vectorBins*rU.refLength(),
                        rhoM_*rU.refMassDensity()
                    );
                }
                    // output velocities
                if(outputField_[1])
                {
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_SAM_SI.xyz",
                        bins*rU.refLength(),
                        USAM_*rU.refVelocity()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_SAM_3D_pos_SI.xyz",
                        vectorBins*rU.refLength(),
                        USAM_*rU.refVelocity()
                    );    
        
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_CAM_SI.xyz",
                        bins*rU.refLength(),
                        UCAM_*rU.refVelocity()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_U_CAM_3D_pos_SI.xyz",
                        vectorBins*rU.refLength(),
                        UCAM_*rU.refVelocity()
                    );   
                }
                // output temperature
                if(outputField_[2])
                {
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_T_SI.xy",
                        bins*rU.refLength(),
                        T_*rU.refTemp()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_T_3D_pos_SI.xy",
                        vectorBins*rU.refLength(),
                        T_*rU.refTemp()
                    );  
                }
                // output pressure
                if(outputField_[3])
                {
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_stress_SI.xyz",
                        bins*rU.refLength(),
                        stress_*rU.refPressure()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_stress_3D_pos_SI.xyz",
                        vectorBins*rU.refLength(),
                        stress_*rU.refPressure()
                    );  
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_p_SI.xy",
                        bins*rU.refLength(),
                        p_*rU.refPressure()
                    );
        
                    writeTimeData
                    (
                        timePath_,
                        "bins_OneDim_"+regionName_+"_"+fieldName_+"_p_3D_pos_SI.xy",
                        vectorBins*rU.refLength(),
                        p_*rU.refPressure()
                    );
                }
            }
        }
    }
}

void polyBinsMethod::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyBinsMethod::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyBinsMethod::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
