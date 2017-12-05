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

#include "polyDensityRadialNew.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyDensityRadialNew, 0);

addToRunTimeSelectionTable(polyField, polyDensityRadialNew, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// void polyDensityRadialNew::setBoundBox
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
polyDensityRadialNew::polyDensityRadialNew
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
//     binWidthX_(lengthX_/scalar(nBinsX_)),
    binWidthY_(lengthY_/scalar(nBinsY_)),
    molIds_(),
    nAvTimeSteps_(0.0),
    resetAtOutput_(true)    
{
//     treshold_ = 1.0;

   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();

    
    // Variable radius
    
    scalar& r=lengthX_;
    scalar Aav = r*r*constant::mathematical::pi/nBinsX_;    
    
    
    DynamicList<scalar> radii;
    DynamicList<scalar> binWidthsX;
    
//     scalar r0=sqrt(Aav/constant::mathematical::pi);
//     radii.append(r0);
//     binWidthsX.append(r0);    

    
    for(label i = 0; i < nBinsX_; i++)
    {
        scalar r1 = 0;
        scalar r2 = 0;
        
        if( i == 0)
        {
            r1=0.0;
            r2=sqrt(Aav/constant::mathematical::pi);
        }
        else
        {
            r1 = radii[i-1];            
            r2 = sqrt((Aav/constant::mathematical::pi) + r1*r1);
        }
        
//         if(r2 <= r)
        {
            radii.append(r2);
            binWidthsX.append(r2-r1);
        }
/*        else
        {
            break;
        }     */   
    }
    
    nBinsX_ = radii.size();
    
    
    
    Info << "Modifying nBinsX to = " << nBinsX_ << endl;
    
    radii_.setSize(nBinsX_);
    binWidthsX_.setSize(nBinsX_);
    
    radii_.transfer(radii);
    binWidthsX_.transfer(binWidthsX);
    
    Info << "radii = " << radii_ << endl;    

    Info << "binWidthsX_ = " << binWidthsX_ << endl;

    r = radii_[nBinsX_-1];
    
    Info << "Modifying radius or lengthX to = " << r << endl;
    
//     Info << "binWidthX = " << binWidthX_ 
//     << ", binWidthY_ = " << binWidthY_
//     << endl;
    
    volumes_.setSize(nBinsX_);
    mass_.setSize(nBinsX_);
    mom_.setSize(nBinsX_);

    rhoM_.setSize(nBinsX_);
    UCAM_.setSize(nBinsX_);
    
    forAll(mass_, i)
    {
        volumes_[i].setSize(nBinsY_, 0.0);
        
        mass_[i].setSize(nBinsY_, 0.0);
        mom_[i].setSize(nBinsY_, vector::zero);
        
        rhoM_[i].setSize(nBinsY_, 0.0);
        UCAM_[i].setSize(nBinsY_, vector::zero);
    }
    
    
    // set volumes and bins
    
    binsX_.setSize(nBinsX_);
    binsY_.setSize(nBinsY_);
    
    forAll(volumes_, i)
    {
        scalar rIn = 0.0;
        scalar rOut = 0.0;
        
        if(i == 0)
        {
            rIn = 0;
            rOut = radii_[i];
        }   
        else
        {
            rIn = radii_[i-1];
            rOut = radii_[i];
        }
        
        binsX_[i] = (rOut+rIn)*0.5; // midpoint of bin
        
        forAll(volumes_[i], j)
        {
            binsY_[j] = binWidthY_*0.5 + binWidthY_*j;
            volumes_[i][j] = constant::mathematical::pi*( (rOut*rOut) - (rIn*rIn) )*binWidthY_;
        }
    }
    
//     Info << "binsX = " << binsX_ << endl;
    
    // framework to make table searching easier
    
    minBinWidth_ = binWidthsX_[nBinsX_-1]/3.0;

    n_.setSize(label(r/minBinWidth_), -1);
    
    forAll(n_, n)
    {
        scalar rI = 0.5*minBinWidth_ + scalar(n)*minBinWidth_;

        for(label i = 0; i < nBinsX_; i++)
        {
            scalar r1 = 0;
            scalar r2 = 0;
            
            if(i == 0)
            {    
                r1 = 0;
                r2 = radii_[i];
            }
            else
            {
                r1 = radii_[i-1];
                r2 = radii_[i];
            }
            
            if((rI >= r1) && (rI < r2))
            {
                n_[n] = i;
            }
        }
    }    
    
//     Info << "n = " << n_ << endl;

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

polyDensityRadialNew::~polyDensityRadialNew()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyDensityRadialNew::createField()
{}

List<label> polyDensityRadialNew::isPointWithinBin
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
    
    label nY = label(rDy/binWidthY_);
    
    if
    (
        (nY >= 0) && (nY < nBinsY_) 
    )
    {
        binNumbers[1] = nY;            
    }    
    
    // radius 
    
    label nX = label(rDR/minBinWidth_);

    if(nX < 0)
    {
        nX = 0;
    }

    if(nX < n_.size())
    {
        scalar i = n_[nX];

        binNumbers[0] = i;
    }
    
    if(binNumbers[0] == -1)
    {
       binNumbers[1] = -1;
    }

    if(binNumbers[1] == -1)
    {
       binNumbers[0] = -1;
    }    
    
    return binNumbers;
}

void polyDensityRadialNew::centreOfMass()
{
//     oldCentreOfMass_

    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    // step 1 - find centre of mass
    vector centre = vector::zero;
    scalar nMols = 0.0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
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

void polyDensityRadialNew::calculateField()
{
    nAvTimeSteps_ += 1.0;
    
    centreOfMass();
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            List<label> n = isPointWithinBin(mol().position());
            
            if((n[0]+n[1]) >= 0)
            {
                const scalar& massI = molCloud_.cP().mass(mol().id());

                mass_[n[0]][n[1]] += massI;
                mom_[n[0]][n[1]] += massI*mol().v();
                
                // measure temperature and so on if required
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
            forAll(mass[i], j)
            {       
                scalar volume = volumes_[i][j];                
                
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


void polyDensityRadialNew::writeToStorage()
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

bool polyDensityRadialNew::readFromStorage()
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

void polyDensityRadialNew::writeField()
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
        }
    }
}

void polyDensityRadialNew::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void polyDensityRadialNew::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& polyDensityRadialNew::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
