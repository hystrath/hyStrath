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

#include "polyFrictionForce.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyFrictionForce, 0);

addToRunTimeSelectionTable(polyField, polyFrictionForce, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// void polyFrictionForce::setBoundBox
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


void polyFrictionForce::setRadius()
{
    // find the radius
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    scalar R = 0.0;
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIdsFluid_, mol().id()) != -1 )
        {
            vector rIS = mol().position() - startPoint_;
            scalar rISD = rIS & n_;
            
            if( (rISD >= 0) && (rISD <= hS_) ) // inside cylinder
            {
                // largest radius?
                scalar rD = mag(rIS - (rISD * n_))
                
                if(rD > R)
                {
                    R = rD;
                }
            }
        }
    }  
    
    // parallel processing
    
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
                    toNeighbour << R;
                }
            }
        }

        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                scalar Rproc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour >> Rproc;
                }
                
                if(Rproc > R)
                {
                    R = Rproc;
                }                
            }
        }
    }    
    
    radius = R + dr_; // plus offset
}

void polyFrictionForce::setCentreOfMass()
{
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
    
    vector rSC = (centre - startPoint_);
    vector rNC = (rSC & n_)* n_;
    vector rSN = rSC - rNC;
    startPoint_ += rSN;
    endPoint_ += rSN;
    endPointS_ += rSN;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyFrictionForce::polyFrictionForce
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
    endPoint_(propsDict_.lookup("endPoint")),
    endPointS_(propsDict_.lookup("endPointSmall")),
    
    
    nMolsMin_(readLabel(propsDict_.lookup("nMolsMin"))),
    dr_(readScalar(propsDict_.lookup("dr"))),
    
    
//     radius_(readScalar(propsDict_.lookup("radius"))),
//     nBins_(readLabel(propsDict_.lookup("nBins"))),
        
//     binWidth_(radius_/scalar(nBins_)),
    h_(mag(endPoint_ - startPoint_)),
    hS_(mag(endPointS_ - startPoint_)),
    n_((endPoint_ - startPoint_)/mag(endPoint_ - startPoint_)),
    

//     magRadii_(),
//     binWidths_(),
//     volumes_(),
//     avVolume_(0.0),
//     minBinWidth_(0.0),
//     n_()    
{
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
    
}




void polyFrictionForce::refreshBins()
{
    setCentreOfMass();
    
    setRadius();
    
    
    
    //- set volumes --- here we impose a constant volume in all bins,
    // hence the binWidth is going change radially

   avVolume_ = radius_*radius_*constant::mathematical::pi*h_/nBins_;

//     Info << "average volume: " << avVolume_ << endl;

    DynamicList<scalar> radii(0);
    DynamicList<scalar> binWidths(0);

    binWidths.append(sqrt(avVolume_/(constant::mathematical::pi*h_)));

    radii.append(0.5*binWidths[0]);

    for(label i = 1; i < nBins_; i++)
    {
        scalar r2 = radii[i-1] + 0.5*binWidths[i-1];
        scalar r1 = sqrt((avVolume_/(constant::mathematical::pi*h_)) + (r2*r2));

        if(r2 <= radius_)
        {
            radii.append(0.5*(r1+r2));
            binWidths.append(r1-r2);
        }
        else
        {
            break;
        }
    }

    //binWidths.shrink();
    //radii.shrink();

    nBins_ = binWidths.size();

    magRadii_.setSize(nBins_, 0.0);
    binWidths_.setSize(nBins_, 0.0);
//     volumes_.setSize(nBins_, 0.0);

    forAll(magRadii_, n)
    {
        magRadii_[n] = radii[n];
        binWidths_[n] = binWidths[n];

//         if(n == 0)
//         {
//             volumes_[n] = mathematicalConstant::pi*binWidths_[n]*binWidths_[n]*h_;
//         }
//         else
//         {
//             volumes_[n] = 2.0*mathematicalConstant::pi*magRadii_[n]*binWidths_[n]*h_;
//         }
    }

//     Info << "volumes: " << volumes_ << " avVolume: " << avVolume_ << endl;

    minBinWidth_ = binWidths_[nBins_-1]/3.0;

    n_.setSize(label(radius_/minBinWidth_), 0);

    forAll(n_, n)
    {
        scalar r = 0.5*minBinWidth_ + scalar(n)*minBinWidth_;

        for(label i = 0; i < nBins_; i++)
        {
            scalar rLimit1 = magRadii_[i] - 0.5*binWidths_[i];
            scalar rLimit2 = magRadii_[i] + 0.5*binWidths_[i];

            if((r >= rLimit1) && (r < rLimit2))
            {
                n_[n] = i;
            }
        }
    }

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyFrictionForce::~polyFrictionForce()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyFrictionForce::createField()
{
    refreshBins();
    
    
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

List<label> polyFrictionForce::isPointWithinBin
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



void polyFrictionForce::calculateField()
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
                    FatalErrorIn("polyFrictionForce::polyFrictionForce()")
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


void polyFrictionForce::writeToStorage()
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

bool polyFrictionForce::readFromStorage()
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

void polyFrictionForce::writeField()
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

void polyFrictionForce::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

    
void polyFrictionForce::measureDuringForceComputationSite
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


const propertyField& polyFrictionForce::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
