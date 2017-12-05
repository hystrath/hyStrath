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

#include "polyDensityFadeBounded.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDensityFadeBounded, 0);
addToRunTimeSelectionTable(polyStateController, polyDensityFadeBounded, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyDensityFadeBounded::polyDensityFadeBounded
(
    Time& t,
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyStateController(t, molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    samplingVolume_(0.0),
    controlVolume_(0.0),
	pointsBox_(molCloud_.rndGen()), 
	molId_(-1),
    maxN_(readLabel(propsDict_.lookup("maxNumberMols"))),
    insertionScheme_(molCloud.rndGen()), //Construct polyFadeII object from polyCloud cachedRandomMD object
	deltaN_(false),
	molsToControlN_(0),
	molsControlled_(0),
	mols_(0.0),
	nMeas_(0.0),
	rhoN_(0.0),
	controlTimeIndex_(0),
	nControlSteps_(0),
	random_(false),
	uniform_(true),
	output_(false)
{

    writeInTimeDir_ = true;
    writeInCase_ = true;

//     singleValueController() = true;

    setBoundBox(propsDict_, bbsampling_, "samplingBoundBox");
    setBoundBox(propsDict_, bbcontrol_, "controlBoundBox");
    samplingVolume_ = bbsampling_.volume();
    controlVolume_ = bbcontrol_.volume();
   
    //- select polyMolecule id for insertion
    {
        const List<word>& idList(molCloud_.cP().molIds());
        const word molId = propsDict_.lookup("molId");
        molId_ = findIndex(idList, molId);
    
        if(molId_ == -1)
        {
            FatalErrorIn("polyDensityFadeBounded::polyDensityFadeBounded()")
                << "Cannot find molId: " << molId << nl << "in: "
                << time_.time().system()/"controllersDict"
                << exit(FatalError);
        }
    }

    pointsBox_.setBoundedBox(bbcontrol_);
    
    scalar writeInterval = readScalar(t.controlDict().lookup("writeInterval"));
    scalar deltaT = t.deltaT().value();
    
    nControlSteps_ = label((writeInterval/deltaT)+0.5);
    
    Info << nl << "NO. OF CONTROL STEPS = " << nControlSteps_ << endl;
    
    if (propsDict_.found("output"))
    {
        output_ = Switch(propsDict_.lookup("output"));
    }
    
    if (propsDict_.found("distribution"))
    {
        const word option = propsDict_.lookup("distribution");
        uniform_ = false;
        
        if(option == "uniform")
        {
            uniform_ = true;
            Info << "Choosing distribution: UNIFORM" << endl;             
        }
        else if(option == "random")
        {
            random_ = true;
            Info << "Choosing distribution: RANDOM" << endl; 
        }        
    }    
    else
    {
        Info << "Choosing default distribution: UNIFORM" << endl; 
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyDensityFadeBounded::~polyDensityFadeBounded()
{}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * *  * * //

void polyDensityFadeBounded::setBoundBox
(
    const dictionary& propsDict,
    boundedBox& bb,
    const word& name 
)
{
    const dictionary& dict(propsDict.subDict(name));
    
    vector startPoint = dict.lookup("startPoint");
    vector endPoint = dict.lookup("endPoint");

    bb.resetBoundedBox(startPoint, endPoint);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyDensityFadeBounded::initialConfiguration()
{
    Info << "polyDensityFadeBounded: initial configuration" << endl;

//     const polyMolecule::constantProperties& constProp = molCloud_.constProps(molId_);
    const scalar& massI = molCloud_.cP().mass(molId_);
    molMass_ = massI;    
    
    velocity_ = propsDict_.lookup("velocity");
    temperature_ = readScalar(propsDict_.lookup("temperature"));

    if(propsDict_.found("numberDensity"))
    {
        density_ = readScalar(propsDict_.lookup("numberDensity"));
    }
    else if (propsDict_.found("massDensitySI"))
    {
        scalar massDensity = readScalar(propsDict_.lookup("massDensitySI"));
        density_ = massDensity/(molCloud_.redUnits().refMassDensity()*molMass_);

        Info<< "Target mass density: " << massDensity << " (SI), " 
             << massDensity/(molCloud_.redUnits().refMassDensity())
             << " (red. units), number density " 
             << density_ << " (red. units)."
             << endl;
    }
    else if (propsDict_.found("deltaN"))
    {
        deltaN_ = true;
        molsToControlN_ = readLabel(propsDict_.lookup("deltaN"));
        molsControlled_=0;
    }
    else
    {
        scalar massDensity = readScalar(propsDict_.lookup("massDensity"));
        density_ = massDensity/molMass_;
    }
    
    // set initial  properties
    
    const word insertionMethod = propsDict_.lookup("insertionMethod");
    const word deletionMethod = propsDict_.lookup("deletionMethod");    
    
    insertionScheme_.createInitialConfiguration
    (
        propsDict_,
        molId_,
        time_.time().deltaT().value(),
        insertionMethod,
        deletionMethod
    );
    
    insertionScheme_.temperature() = temperature_;
    insertionScheme_.velocity() = velocity_; 
    
    insertionScheme_.output(time_);
    
    if (propsDict_.found("controlFromStart"))
    {
        bool controlFromStart = false;
        
        controlFromStart = Switch(propsDict_.lookup("controlFromStart"));
        
        if(controlFromStart)
        {
            controlTimeIndex_ = nControlSteps_ + 1;
        }
    }
        
}

void polyDensityFadeBounded::controlBeforeVelocityI()
{}

void polyDensityFadeBounded::controlBeforeMove()
{}

void polyDensityFadeBounded::controlBeforeForces()
{}

void polyDensityFadeBounded::controlDuringForces
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}



void polyDensityFadeBounded::controlAfterForces()
{
    insertionScheme_.checkFractions(molCloud_);
    
    controlTimeIndex_++;

    // sampling
    measureProperties();
    
    // control on the first time-step of the write interval
    if(controlTimeIndex_ > nControlSteps_)
    {
        Info << "polyDensityFadeBounded: control ... " << endl;
        
        averageProperties();
        
        label molsToControl = noOfMolsToControl();
        
        label totalMols = molsToControl;
        
        if(mag(totalMols) > 0)
        {        
            nMolsToInsert(molsToControl);
            
            Pout << "mols to control = " << molsToControl << endl; 
            
            List<vector> molPositions(molPoints_.size(), vector::zero);
        
            forAll(molPoints_, i)
            {
                molPositions[i] = molPoints_[i];
            }
        
            insertionScheme_.controlMolecules
            (
                molCloud_,
                molsToControl,
                bbcontrol_,
                molPositions 
            );

            if(deltaN_)
            {
                label nInserted = insertionScheme_.nInserted();
                label nDeleted = insertionScheme_.nDeleted();
        
                molsControlled_ += nInserted;
                molsControlled_ -= nDeleted;
            } 
        }

        mols_ = 0.0;
        nMeas_ = 0.0;
        controlTimeIndex_ = 1;
    }
}

void polyDensityFadeBounded::measureProperties()
{
    scalar mols = 0.0;
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());

    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(bbsampling_.contains(mol().position()))
        {
            const label& tN = mol().trackingNumber();

            if(mol().id() == molId_)
            {
                if(findIndex(insertionScheme_.deletingList(), tN) == -1)
                {
                    mols += 1.0;
                }
            }
        }
    }
    
    if(Pstream::parRun())
    {
        reduce(mols, sumOp<scalar>());
    }

    mols_ += mols;
    nMeas_ += 1.0;
}

void polyDensityFadeBounded::averageProperties()
{
    rhoN_ = mols_/(nMeas_*samplingVolume_);
    
    if(Pstream::parRun())
    {
        Pout<< "[measured] number density (R.U.) = " << rhoN_ 
            << ", mass density (R.U.) = " << mols_*molMass_/(nMeas_*samplingVolume_)
            << ", mass density (S.I.) = " 
            << molCloud_.redUnits().refMassDensity()*mols_*molMass_/(nMeas_*samplingVolume_)
            << endl;
    }
}


label polyDensityFadeBounded::noOfMolsToControl()
{
    scalar nMols = (density_ - rhoN_)*samplingVolume_;

    if(deltaN_)
    {
        nMols = molsToControlN_ - molsControlled_;
        
        Info << "[Delta N]"
             << " Target mols to control = " << molsToControlN_ 
             << ", mols controlled = " << molsControlled_
             << ", controlling left = " <<  nMols
             << endl;
             
    }

    if(mag(nMols) > maxN_)
    {
        if(nMols > 0.0)
        {
            nMols = maxN_;
        }
        else if(nMols < 0.0)
        {
            nMols = -maxN_;
        }
    }

    // -total number of polyMolecules to insert within entire control zone
    label nMolsRounded = 0;

    if(nMols > 0.0)
    {
        nMolsRounded = label(nMols + 0.5);
    }
    else if(nMols < 0.0)
    {
        nMolsRounded = label(nMols - 0.5);
    }

    return nMolsRounded;
}

void polyDensityFadeBounded::nMolsToInsert(label& molsToControl)
{
    vectorField nMolsPoints(mag(molsToControl), vector::zero);

    if(Pstream::master())
    {
        if(random_)
        {
            //randomly generate points
            forAll(nMolsPoints, p)
            {
                nMolsPoints[p] = pointsBox_.randomPoint();
            }
        }
        else if(uniform_)
        {
            nMolsPoints = pointsBox_.uniform(mag(molsToControl));
        }        
    }
        
    if(Pstream::parRun())
    {
        forAll(nMolsPoints, i)
        {
            reduce(nMolsPoints[i], sumOp<vector>());         
        }
    }
 
    DynamicList<vector> molPoints;     
    
    List<bool> pointsOnThisProc(mag(molsToControl), false);
    
    forAll(nMolsPoints, i)
    {
        label cellI = mesh_.findCell(nMolsPoints[i]);

        if(cellI != -1)
        {
            molPoints.append(nMolsPoints[i]);
            pointsOnThisProc[i]=true;
        }
    } 
 
    // test for bug - are same molecules on different processors
    if(Pstream::parRun())
    { 
        //- sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);
                    toNeighbour << pointsOnThisProc;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<bool> pointsOnThisProcRec;
                
                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> pointsOnThisProcRec;
                }
                
                forAll(pointsOnThisProcRec, i)
                {
                    if(pointsOnThisProcRec[i] && pointsOnThisProc[i])
                    {
                        if(Pstream::myProcNo() > p)
                        {
                            pointsOnThisProc[i]=false;
                        }    
                    }
                }
            }
        }
        
        molPoints.clear();
        
        forAll(nMolsPoints, i)
        {
            if(pointsOnThisProc[i])
            {
                molPoints.append(nMolsPoints[i]);
            }
        }
    }
    
    //molPoints.shrink();
    molPoints_.clear();
    molPoints_.transfer(molPoints);
    
    label newMolsToControl = 0;    
    
    if(molsToControl > 0)
    {
        newMolsToControl = molPoints_.size();
    }
    else if (molsToControl < 0)
    {
        newMolsToControl = -molPoints_.size();
    }
    
    if(Pstream::parRun())
    {
        label totalMols = newMolsToControl;
          
        reduce(totalMols, sumOp<label>());

        if(mag(molsToControl) > 0)
        {
            Pout<< "measured number density: " << rhoN_ 
                <<", total nMols (target): " << molsToControl 
                << ", total nMols (check): " << totalMols
                << ", this processor nMols: "<< newMolsToControl
                << endl;
        }
    }
    else
    {
        label totalMols = newMolsToControl;

        if(mag(molsToControl) > 0)
        {
            Info<< "measured density: " << rhoN_ 
                <<", total nMols (target): " << molsToControl 
                << ", total nMols (check): " << totalMols
                << endl;
        }
    }
    
    molsToControl = newMolsToControl;    
}


void polyDensityFadeBounded::controlAfterVelocityII()
{

}

void polyDensityFadeBounded::calculateProperties()
{}

void polyDensityFadeBounded::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}

void polyDensityFadeBounded::updateProperties(const dictionary& newDict)
{
    //- the main controller properties should be updated first
    updateStateControllerProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    insertionScheme_.updateProperties(propsDict_);
}


} // End namespace Foam

// ************************************************************************* //
