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

#include "polyForceZone.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyForceZone, 0);
addToRunTimeSelectionTable(polyField, polyForceZone, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyForceZone::setBoundBoxes()
{
    PtrList<entry> boxList(propsDict_.lookup("boxes"));

    boxes_.setSize(boxList.size());

    forAll(boxList, b)
    {
        const entry& boxI = boxList[b];
        const dictionary& dict = boxI.dict();

        vector startPoint = dict.lookup("startPoint");
        vector endPoint = dict.lookup("endPoint");
        boxes_[b].resetBoundedBox(startPoint, endPoint);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyForceZone::polyForceZone
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
    boxes_(),
    totalVolume_(0.0),
    force_(vector::zero),
    molIds_(),
    timeIndex_(0),
	forceField_(1, vector::zero),

    nAvTimeSteps_(0.0),
    resetAtOutput_(true)    
{
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
        nameFile_ = "propertiesZoneBoundedData_"+fieldName_;

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
    
    // build bound boxes

    setBoundBoxes();

    //-set the total volume
    forAll(boxes_, b)
    {
        totalVolume_ += boxes_[b].volume();
    }

    // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyForceZone::~polyForceZone()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyForceZone::createField()
{}

void polyForceZone::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
//     vector force = vector::zero;
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());

        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIds_, mol().id()) != -1)
            {
                forAll(boxes_, b)
                {
                    if(boxes_[b].contains(mol().position()))
                    {
//                         const scalar& massI = molCloud_.cP().mass(mol().id());

                        forAll(mol().siteForces(), i)
                        {
                            force_ += mol().siteForces()[i];
                            
//                             Pout << "force = " << mol().siteForces()[i] << endl;
                        }
                    }
                }
            }
        }
    }

//     force_ += force;

    if(time_.outputTime())
    {
        //- parallel communication
        vector force = force_;

        if(Pstream::parRun())
        {
            reduce(force, sumOp<vector>());
        }

        const scalar& nAvTimeSteps = nAvTimeSteps_;

        // density

        forceField_[timeIndex_] = force/nAvTimeSteps;

        if(resetAtOutput_)
        {
            nAvTimeSteps_ = 0.0;
            force_ = vector::zero;
        }
        else 
        {
            writeToStorage();
        }
    }
}

void polyForceZone::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
        file << nAvTimeSteps_ << endl;
        file << force_ << endl;
    }
    else
    {
        FatalErrorIn("void polyForceZone::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool polyForceZone::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar nAvTimeSteps;
        vector force;

        file >> nAvTimeSteps;
        file >> force;

        nAvTimeSteps_ = nAvTimeSteps;
        force_ = force;
    }

    return goodFile;    
}

void polyForceZone::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeField(1, runTime.timeOutputValue());
            
            writeTimeData
            (
                casePath_,
                "forceZone_bb_"+fieldName_+"_F.xyz",
                timeField,
                forceField_,
                true
            );

            const reducedUnits& rU = molCloud_.redUnits();
    
            if(rU.outputSIUnits())
            {
                writeTimeData
                (
                    casePath_,
                    "forceZone_bb_"+fieldName_+"_F_SI.xyz",
                    timeField*rU.refTime(),
                    forceField_*rU.refForce(),
                    true
                );

            }
        }
    }
}


void polyForceZone::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void polyForceZone::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& polyForceZone::fields() const
{
    return fields_;
}

} // End namespace Foam

// ************************************************************************* //
