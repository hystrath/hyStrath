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

#include "poly2DBinsSites.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(poly2DBinsSites, 0);
addToRunTimeSelectionTable(polyField, poly2DBinsSites, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
poly2DBinsSites::poly2DBinsSites
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

    molSiteIds_(),

    atoms_(),

    N_(),
    

    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("poly2DBinsSites::poly2DBinsSites()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
    
    // choose molecule ids to sample

    molSiteIds_.clear();

    selectSiteIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molSiteIds_ = ids.siteIds();


    
    // create bin model
    binModel_ = autoPtr<twoDimBinModel>
    (
        twoDimBinModel::New(mesh, propsDict_)
    );
    
    List<label> nBins = binModel_->nBins();
    
    label nBinsX = nBins[0];
    label nBinsY = nBins[1];

//     Info << "nBinsX = " << nBinsX << endl;
//     Info << "nBinsY = " << nBinsY << endl;
    
    atoms_.setSize(nBinsX);

    N_.setSize(nBinsX);
    
    forAll(atoms_, i)
    {
        atoms_[i].setSize(nBinsY, 0.0);
        N_[i].setSize(nBinsY, 0.0);
    }
    
//     Info << "atoms_ = " << atoms_ << endl;
//     Info << "mom_ = " << mom_ << endl;        
    
   
    
    // read in stored data from dictionary

    resetAtOutput_ = Switch(propsDict_.lookup("resetAtOutput"));
    
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
   
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

poly2DBinsSites::~poly2DBinsSites()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void poly2DBinsSites::createField()
{
}


void poly2DBinsSites::calculateField()
{
//     Info << "atoms = " << atoms_ << endl;
//     Info << "mom = " << mom_ << endl;
    
    nAvTimeSteps_ += 1.0;
    
    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];
            
            forAll(molI->sitePositions(), i)
            {
                label siteId = molCloud_.cP().siteNames_to_siteIdList()[molI->id()][i];
        
                if(findIndex(molSiteIds_, siteId) != -1)
                {
//                     Info<< "siteId = " << siteId 
//                         << ", in siteIdList = " << molCloud_.cP().siteIds()[siteId]
//                         << endl;
                    
                    const vector& rI = molI->sitePositions()[i];
                    
                    List<label> n = binModel_->isPointWithinBin(rI, cellI);

                    if((n[0]+n[1]) >= 0)
                    {
                        atoms_[n[0]][n[1]] += 1.0;
                    }
                }
            }

        }
    }

    if(time_.outputTime()) 
    {
        List<scalarField> atoms = atoms_;
       

        if(Pstream::parRun())
        {
            forAll(atoms, i)
            {
                forAll(atoms[i], j)
                {
                    reduce(atoms[i][j], sumOp<scalar>());
                }
            }
        }
        
        const scalar& nAvTimeSteps = nAvTimeSteps_;
        
        forAll(atoms, i)
        {
//             scalar volume = binModel_->binVolume(i);
            
            forAll(atoms[i], j)
            {       
                N_[i][j] = atoms[i][j]/(nAvTimeSteps);

            }
        }
    


        if(resetAtOutput_)
        {
            //- reset fields
            nAvTimeSteps_ = 0.0;
            
            forAll(atoms, i)
            {
                forAll(atoms[i], j)
                {
                    atoms_[i][j]=0.0;
                }
            }
        }
        else //if(averagingAcrossManyRuns_)
        {
            writeToStorage();
        }
    }
}

void poly2DBinsSites::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
        file << nAvTimeSteps_ << endl;
        file << atoms_ << endl;
//         file << mom_ << endl; 
    }
    else
    {
        FatalErrorIn("void poly2DBinsSites::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool poly2DBinsSites::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar nAvTimeSteps;

        file >> nAvTimeSteps;
        file >> atoms_;
//         file >> mom_;

        nAvTimeSteps_ = nAvTimeSteps;
    }

    return goodFile;
}

void poly2DBinsSites::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            
            scalarField binsX = binModel_->binPositionsX();
            scalarField binsY = binModel_->binPositionsY();            
            

            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_binsX.xy",
                binsX
            ); 
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_binsY.xy",
                binsY
            );           

            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_N.xy",
                N_
            );           
            
            
            
            const reducedUnits& rU = molCloud_.redUnits();
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_binsX_SI.xy",
                binsX*rU.refLength()
            ); 
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_binsY_SI.xy",
                binsY*rU.refLength()
            );           
        }
    }
}

void poly2DBinsSites::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void poly2DBinsSites::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& poly2DBinsSites::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
