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

#include "polyBindingEnergyRadial.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyBindingEnergyRadial, 0);

addToRunTimeSelectionTable(polyField, polyBindingEnergyRadial, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// void polyBindingEnergyRadial::setBoundBox
// (
//     const dictionary& propsDict,
//     boundedBox& bb,
//     const word& name 
// )
// {
//     const dictionary& dict(propsDict.subDict(name));
//     
//     vector startPoint = dict.lookup("startPoint");
//     vector endPoint = dict.lookup("endPoint");
//     bb.resetBoundedBox(startPoint, endPoint);
// }




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyBindingEnergyRadial::polyBindingEnergyRadial
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
    
    startPoint_(propsDict_.lookup("startPoint")),
    unitVectorX_(propsDict_.lookup("unitVectorX")),
    unitVectorY_(propsDict_.lookup("unitVectorY")),
    unitVectorZ_(propsDict_.lookup("unitVectorZ")),
    nBinsX_(readLabel(propsDict_.lookup("nBinsX"))),
    nBinsY_(readLabel(propsDict_.lookup("nBinsY"))),
    lengthX_(readLabel(propsDict_.lookup("lengthX"))),    
    lengthY_(readLabel(propsDict_.lookup("lengthY"))),
    binWidthX_(lengthX_/scalar(nBinsX_)),
    binWidthY_(lengthY_/scalar(nBinsY_)),
    resetAtOutput_(true)    
{
//     treshold_ = 1.0;
    measureInterForcesSites_ = true;
    
    // choose molecule ids to sample
    {
        // choose molecule ids to sample
        molIdsWall_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsWall"
            
        );

        molIdsWall_ = ids.molIds();
    }
    
    {
        // choose molecule ids to sample
        molIdsFluid_.clear();
        
        selectIds ids
        (
            molCloud_.cP(),
            propsDict_,
            "molIdsFluid"
            
        );

        molIdsFluid_ = ids.molIds();
    }

    Info << "binWidthX = " << binWidthX_ 
    << ", binWidthY_ = " << binWidthY_
    << endl;
    
    volumes_.setSize(nBinsX_);
    
    pairs_.setSize(nBinsX_);
    mols_.setSize(nBinsX_);
    energy_.setSize(nBinsX_);
    
    bindingEnergy_.setSize(nBinsX_);
    avNoOfPairs_.setSize(nBinsX_);
    
    forAll(mols_, i)
    {
        volumes_[i].setSize(nBinsY_, 0.0);
        
        pairs_[i].setSize(nBinsY_, 0.0);
        mols_[i].setSize(nBinsY_, 0.0);
        energy_[i].setSize(nBinsY_, 0.0);
        
        bindingEnergy_[i].setSize(nBinsY_, 0.0);
        avNoOfPairs_[i].setSize(nBinsY_, 0.0);
    }
    
    
    // set volumes and bins
    
    binsX_.setSize(nBinsX_);
    binsY_.setSize(nBinsY_);
    
    forAll(volumes_, i)
    {
        scalar rIn = binWidthX_*i;
        scalar rOut = binWidthX_*(i+1);
        
        binsX_[i] = binWidthX_*0.5 + binWidthX_*i;
        
        forAll(volumes_[i], j)
        {
            binsY_[j] = binWidthY_*0.5 + binWidthY_*j;
            volumes_[i][j] = constant::mathematical::pi*( (rOut*rOut) - (rIn*rIn) )*binWidthY_;
        }
    }

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
//             Info << "Reading from storage, e.g. noAvTimeSteps = " << nAvTimeSteps_ << endl;
        }
    }    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyBindingEnergyRadial::~polyBindingEnergyRadial()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyBindingEnergyRadial::createField()
{
    
    label N = molCloud_.moleculeTracking().getMaxTrackingNumber();
    
    Pout << "N = " << N << endl;
    
    molPositions_.setSize(N, vector::zero);
    fluidMols_.setSize(N, false);
    nPairs_.setSize(N, 0.0);
    DeltaE_.setSize(N, 0.0);
    
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIdsFluid_, mol().id()) != -1 )
        {
            label tN = mol().trackingNumber();
            fluidMols_[tN] = true;
        }
    }     
    
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
                    toNeighbour << fluidMols_;
                }
            }
        }

        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {

                List<bool> fluidMolsProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> fluidMolsProc;
                }
                
                forAll(fluidMolsProc, i)
                {
                    if(fluidMolsProc[i])
                    {
                        fluidMols_[i] = true;
                    }
                }
            }
        }
    }
    
}

List<label> polyBindingEnergyRadial::isPointWithinBin
(
    const vector& rI
)
{
    List<label> binNumbers;
binNumbers.append(-1);
binNumbers.append(-1);
    
    vector rSI = rI - startPoint_;
    
    vector unitVectorR = (rSI & unitVectorX_)*unitVectorX_ + (rSI & unitVectorZ_)*unitVectorZ_;

    unitVectorR /= mag(unitVectorR);
    
    scalar rDR = rSI & unitVectorR;
    
    scalar rDy = rSI & unitVectorY_;
    label nX = label(rDR/binWidthX_);
    label nY = label(rDy/binWidthY_);
    
    if
    (
        ( (nX >= 0) && (nY >= 0) ) &&
        ( (nX < nBinsX_) && (nY < nBinsY_) )
    )
    {
        binNumbers[0] = nX;
        binNumbers[1] = nY;            
    }    
    
    return binNumbers;
}

void polyBindingEnergyRadial::centreOfMass()
{
//     oldCentreOfMass_

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    // step 1 - find centre of mass
    vector centre = vector::zero;
    scalar nMols = 0.0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIdsFluid_, mol().id()) != -1)
        {
            centre += mol().position();
            nMols += 1.0;
        }
    }    
    
    if(Pstream::parRun())
    {
        reduce(centre, sumOp<vector>());
        reduce(nMols, sumOp<scalar>());
    }

    if(nMols > 0)
    {
        centre /= nMols;
    }
    
    centreOfMass_ = centre;

    centreOfMass_ = (centreOfMass_ & unitVectorX_) * unitVectorX_
    + (startPoint_ & unitVectorY_)* unitVectorY_
    + (centreOfMass_ & unitVectorZ_)*unitVectorZ_;    
}

void polyBindingEnergyRadial::calculateField()
{
//     nAvTimeSteps_ += 1.0;
    
    centreOfMass();
    
    {
        IDLList<polyMolecule>::iterator mol(molCloud_.begin());
        
        for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
        {
            if(findIndex(molIdsFluid_, mol().id()) != -1 )
            {
                label tN = mol().trackingNumber();
                
                if(tN >= molPositions_.size())
                {
                    FatalErrorIn("polyBindingEnergyRadial::polyBindingEnergyRadial()")
                        << "Oops something went wrong."
                        << exit(FatalError);                
                }
                
                molPositions_[tN] = mol().position();
            }
        }    
    }
    
    if(Pstream::parRun())
    {
        forAll(DeltaE_, i)
        {
            reduce(molPositions_[i], sumOp<vector>());            
            reduce(DeltaE_[i], sumOp<scalar>());
            reduce(nPairs_[i], sumOp<scalar>());
        }
    }    
    
//     Pout << "Fluid Mols = " << fluidMols_ << endl;

    
    forAll(molPositions_, i)
    {
        if(fluidMols_[i])
        {
//             Info << "molPositions_[i] = " << molPositions_[i] << endl;
            
            List<label> n = isPointWithinBin(molPositions_[i]);
            
            if((n[0]+n[1]) >= 0)
            {
                energy_[n[0]][n[1]] += DeltaE_[i];
                mols_[n[0]][n[1]] += 1.0;
                pairs_[n[0]][n[1]] += nPairs_[i];
            }
            
         
        }
        
        // reset 
        DeltaE_[i] = 0.0;
        nPairs_[i] = 0.0;         
        molPositions_[i] = vector::zero;        
    }


    if(time_.outputTime()) 
    {
        forAll(bindingEnergy_, i)
        {
            forAll(bindingEnergy_[i], j)
            {       
                if(mols_[i][j] > 0)
                {
                    bindingEnergy_[i][j] = energy_[i][j]/mols_[i][j];
                    avNoOfPairs_[i][j] = pairs_[i][j]/mols_[i][j];
                }
            }
        }

        if(resetAtOutput_)
        {
            forAll(energy_, i)
            {
                forAll(energy_[i], j)
                {
                    energy_[i][j]=0.0;
                    pairs_[i][j]=0.0;
                    mols_[i][j]=0.0;
                }
            }
        }
        else //if(averagingAcrossManyRuns_)
        {
            writeToStorage();
        }
    }    

}


void polyBindingEnergyRadial::writeToStorage()
{
    OFstream file(pathName_/nameFile_);

    if(file.good())
    {
//         file << nAvTimeSteps_ << endl;
        file << pairs_ << endl;
        file << mols_ << endl;
        file << energy_ << endl; 
    }
    else
    {
        FatalErrorIn("void poly2DBins::writeToStorage()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

bool polyBindingEnergyRadial::readFromStorage()
{
    IFstream file(pathName_/nameFile_);

    bool goodFile = file.good();

    if(goodFile)
    {
        file >> pairs_;
        file >> mols_;
        file >> energy_;
    }

    return goodFile;
}

void polyBindingEnergyRadial::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            const scalarField& binsX = binsX_;
            const scalarField& binsY = binsY_;
            
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
                "bins_twoDim_"+fieldName_+"_bindingEnergy.xy",
                bindingEnergy_
            );     
            
            writeTimeData
            (
                timePath_,
                "bins_twoDim_"+fieldName_+"_noOfPairs.xy",
                avNoOfPairs_
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
                "bins_twoDim_"+fieldName_+"_bindingEnergy_SI.xy",
                bindingEnergy_*rU.refEnergy()
            );              
        }
    }
}

void polyBindingEnergyRadial::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

    
void polyBindingEnergyRadial::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
)
{
    label idI = molI->id();
    label idJ = molJ->id();
    
    label idIF = findIndex(molIdsFluid_, molI->id());
    label idJF = findIndex(molIdsFluid_, molJ->id());

    label idIW = findIndex(molIdsWall_, molI->id());
    label idJW = findIndex(molIdsWall_, molJ->id());

    label tN = -1;
    
    if
    (
        ((idIW != -1) && (idJF != -1)) 
    )
    {
        tN = molJ->trackingNumber();
    }
    else if((idJW != -1) && (idIF != -1))
    {
        tN = molI->trackingNumber();        
    }    
    
    if(tN != -1)
    {
        label k = molCloud_.pot().pairPots().pairPotentialIndex(idI, idJ, sI, sJ);    

        vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
        scalar rsIsJMag = mag(rsIsJ);
        scalar pE = molCloud_.pot().pairPots().energy(k, rsIsJMag);
        
        if(molI->referred() || molJ->referred())
        {
            nPairs_[tN] += 0.5;
            DeltaE_[tN] += pE*0.5*0.5; 
        }
        else
        {
            nPairs_[tN] += 1;
            DeltaE_[tN] += pE*0.5; 
        }
    }
    
}


const propertyField& polyBindingEnergyRadial::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
