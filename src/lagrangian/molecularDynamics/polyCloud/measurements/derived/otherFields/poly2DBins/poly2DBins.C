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

#include "poly2DBins.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(poly2DBins, 0);
addToRunTimeSelectionTable(polyField, poly2DBins, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
poly2DBins::poly2DBins
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

    mass_(),
    mom_(),

    rhoM_(),
    UCAM_(),

    nAvTimeSteps_(0.0),
    resetAtOutput_(true)
{
    
    const cellZoneMesh& cellZones = mesh_.cellZones();
    
    regionId_ = cellZones.findZoneID(regionName_);

    if(regionId_ == -1)
    {
        FatalErrorIn("poly2DBins::poly2DBins()")
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
    binModel_ = autoPtr<twoDimBinModel>
    (
        twoDimBinModel::New(mesh, propsDict_)
    );
    
    List<label> nBins = binModel_->nBins();
    
    label nBinsX = nBins[0];
    label nBinsY = nBins[1];

//     Info << "nBinsX = " << nBinsX << endl;
//     Info << "nBinsY = " << nBinsY << endl;
    
    mass_.setSize(nBinsX);
    mom_.setSize(nBinsX);

    rhoM_.setSize(nBinsX);
    UCAM_.setSize(nBinsX);
    
    forAll(mass_, i)
    {
        mass_[i].setSize(nBinsY, 0.0);
        mom_[i].setSize(nBinsY, vector::zero);
        
        rhoM_[i].setSize(nBinsY, 0.0);
        UCAM_[i].setSize(nBinsY, vector::zero);
    }
    
//     Info << "mass_ = " << mass_ << endl;
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

poly2DBins::~poly2DBins()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void poly2DBins::createField()
{
}


void poly2DBins::calculateField()
{
//     Info << "mass = " << mass_ << endl;
//     Info << "mom = " << mom_ << endl;
    
    nAvTimeSteps_ += 1.0;
    
    forAll(mesh_.cellZones()[regionId_], c)
    {
        const label& cellI = mesh_.cellZones()[regionId_][c];
        const List<polyMolecule*>& molsInCell = molCloud_.cellOccupancy()[cellI];

        forAll(molsInCell, mIC)
        {
            polyMolecule* molI = molsInCell[mIC];

            const vector& rI = molI->position();

            List<label> n = binModel_->isPointWithinBin(rI, cellI);
//             Info << "n = "<< n << endl;
            if((n[0]+n[1]) >= 0)
            {
                if(findIndex(molIds_, molI->id()) != -1)
                {
                    const scalar& massI = molCloud_.cP().mass(molI->id());

                    mass_[n[0]][n[1]] += massI;
                    mom_[n[0]][n[1]] += massI*molI->v();
                }
            }
        }
    }

    if(time_.outputTime()) 
    {
        List<scalarField> mass = mass_;
        List<vectorField> mom = mom_;        

        if(Pstream::parRun())
        {
            forAll(mass, i)
            {
                forAll(mass[i], j)
                {
                    reduce(mass[i][j], sumOp<scalar>());
                    reduce(mom[i][j], sumOp<vector>());
                }
            }
        }
        
        const scalar& nAvTimeSteps = nAvTimeSteps_;
        
        forAll(mass, i)
        {
            scalar volume = binModel_->binVolume(i);
            
            forAll(mass[i], j)
            {       
                rhoM_[i][j] = mass[i][j]/(nAvTimeSteps*volume);
                
                if(mass[i][j] > 0)
                {
                    UCAM_[i][j] = mom[i][j]/mass[i][j];
                }
            }
        }
    


        if(resetAtOutput_)
        {
            //- reset fields
            nAvTimeSteps_ = 0.0;
            
            forAll(mass, i)
            {
                forAll(mass[i], j)
                {
                    mass_[i][j]=0.0;
                    mom_[i][j]=vector::zero;
                }
            }
        }
        else //if(averagingAcrossManyRuns_)
        {
            writeToStorage();
        }
    }
}

void poly2DBins::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
        file << nAvTimeSteps_ << endl;
        file << mass_ << endl;
        file << mom_ << endl; 
    }
    else
    {
        FatalErrorIn("void poly2DBins::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool poly2DBins::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        scalar nAvTimeSteps;

        file >> nAvTimeSteps;
        file >> mass_;
        file >> mom_;

        nAvTimeSteps_ = nAvTimeSteps;
    }

    return goodFile;
}

void poly2DBins::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            
            scalarField binsX = binModel_->binPositionsX();
            scalarField binsY = binModel_->binPositionsY();            
            

            // new write out 
            {
                OFstream file(timePath_/"bins_twoDim_"+fieldName_+"_rhoM_II.xy");

                if(file.good())
                {
                    forAll(binsY, j)
                    {
                        forAll(binsX, i)
                        {
                            file 
                                << binsX[i] << "\t" 
                                << binsY[j] << "\t" 
                                << rhoM_[i][j] << "\t"
                                << endl;
                        }
                    }
                }
                else
                {
                    FatalErrorIn("void writeTimeData::writeTimeData()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }            
            }
            {
                OFstream file(timePath_/"bins_twoDim_"+fieldName_+"_U_CAM_II_X.xy");

                if(file.good())
                {
                    forAll(binsY, j)
                    {
                        forAll(binsX, i)
                        {
                            file 
                                << binsX[i] << "\t" 
                                << binsY[j] << "\t" 
                                << UCAM_[i][j].x() << "\t"
                                << endl;
                        }
                    }
                }
                else
                {
                    FatalErrorIn("void writeTimeData::writeTimeData()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }            
            }      
            
            {
                OFstream file(timePath_/"bins_twoDim_"+fieldName_+"_U_CAM_II_Y.xy");

                if(file.good())
                {
                    forAll(binsY, j)
                    {
                        forAll(binsX, i)
                        {
                            file 
                                << binsX[i] << "\t" 
                                << binsY[j] << "\t" 
                                << UCAM_[i][j].y() << "\t"
                                << endl;
                        }
                    }
                }
                else
                {
                    FatalErrorIn("void writeTimeData::writeTimeData()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }            
            }
            {
                OFstream file(timePath_/"bins_twoDim_"+fieldName_+"_U_CAM_II_Z.xy");

                if(file.good())
                {
                    forAll(binsY, j)
                    {
                        forAll(binsX, i)
                        {
                            file 
                                << binsX[i] << "\t" 
                                << binsY[j] << "\t" 
                                << UCAM_[i][j].z() << "\t"
                                << endl;
                        }
                    }
                }
                else
                {
                    FatalErrorIn("void writeTimeData::writeTimeData()")
                        << "Cannot open file " << file.name()
                        << abort(FatalError);
                }            
            }            
            
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
                "bins_twoDim_"+fieldName_+"_rhoM.xy",
                rhoM_
            );           
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_X.xy",
                UCAM_,
                "x"
            );  
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_Y.xy",
                UCAM_,
                "y"
            ); 

            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_Z.xy",
                UCAM_,
                "z"
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

            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_rhoM_SI.xy",
                rhoM_*rU.refMassDensity()
            );           
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_X_SI.xy",
                UCAM_*rU.refVelocity(),
                "x"
            );  
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_Y_SI.xy",
                UCAM_*rU.refVelocity(),
                "y"
            ); 

            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_U_CAM_Z_SI.xy",
                UCAM_*rU.refVelocity(),
                "z"
            );             

//             scalarField bins = binModel_->binPositionsX();
//             scalarField binsY = binModel_->binPositionsY();
// 
//             forAll(rhoM_, i)
//             {
//                 std::string s;
//                 std::stringstream out;
//                 out << bins[i];
//                 s = out.str();
// 
//                 writeTimeData
//                 (
//                     timePath_,
//                     "bins_twoDim_"+fieldName_+"_"+s+"_rhoM.xy",
//                     binsY,
//                     rhoM_[i]
//                 );
// 
//                 writeTimeData
//                 (
//                     timePath_,
//                     "bins_twoDim_"+fieldName_+"_"+s+"_U_CAM.xyz",
//                     binsY,
//                     UCAM_[i]
//                 );
//                 
//                 const reducedUnits& rU = molCloud_.redUnits();
//     
//                 if(rU.outputSIUnits())
//                 {
//                     writeTimeData
//                     (
//                         timePath_,
//                         "bins_twoDim_"+fieldName_+"_"+s+"_rhoM_SI.xy",
//                         binsY*rU.refLength(),
//                         rhoM_[i]*rU.refMassDensity()
//                     );
// 
//                     writeTimeData
//                     (
//                         timePath_,
//                         "bins_twoDim_"+fieldName_+"_"+s+"_U_CAM_SI.xyz",
//                         binsY*rU.refLength(),
//                         UCAM_[i]*rU.refVelocity()
//                     );
//                 }
//             }
        }
    }
}

void poly2DBins::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
)
{}

void poly2DBins::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{}

const propertyField& poly2DBins::fields() const
{
    return  fields_;
}

} // End namespace Foam

// ************************************************************************* //
