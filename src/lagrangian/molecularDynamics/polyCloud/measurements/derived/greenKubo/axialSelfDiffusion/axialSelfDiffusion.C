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

#include "axialSelfDiffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "graph.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(axialSelfDiffusion, 0);

addToRunTimeSelectionTable(polyField, axialSelfDiffusion, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //






// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
axialSelfDiffusion::axialSelfDiffusion
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
    unitVector_(propsDict_.lookup("unitVector")),    
    nSteps_(readLabel(propsDict_.lookup("nSteps")))  

{
   // choose molecule ids to sample

    molIds_.clear();

    selectIds ids
    (
        molCloud_.cP(),
        propsDict_
    );

    molIds_ = ids.molIds();
    
    nS_ = 0;
    nBatch_ = 0.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

axialSelfDiffusion::~axialSelfDiffusion()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void axialSelfDiffusion::createField()
{
    label maxTN = 0; 
    
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    DynamicList<label> tNs(0);
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        label tN = mol().trackingNumber();
        
        if(findIndex(molIds_, mol().id()) != -1)
        {
            tNs.append(tN);
        }
        
        if(tN > maxTN)
        {
            maxTN = tN;
        }
    }
    
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

                    toNeighbour << maxTN;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                label maxTNProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour  >> maxTNProc;
                }
                
                if(maxTNProc > maxTN)
                {
                    maxTN = maxTNProc;
                }
            }
        }
    }    
   
    Pout << "maxTN = " << maxTN << endl;
    
    //- parallel communication
    if(Pstream::parRun())
    {
        List<label> tNsTransfer(tNs.size());
        
        forAll(tNs, i)
        {
            tNsTransfer[i]=tNs[i];
        }
        
        //-sending
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                const int proc = p;
                {
                    OPstream toNeighbour(Pstream::blocking, proc);

                    toNeighbour << tNsTransfer;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<label> tNsProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour  >> tNsProc;
                }
    
                forAll(tNsProc, i)
                {
                    tNs.append(tNsProc[i]);
                }
            }
        }
    }
    
//     Pout << "error1" << endl;
    
    tNaddress_.setSize(maxTN+1, -1);
    
    if(Pstream::myProcNo() == 0)
    {
        forAll(tNs, i)
        {        
            tNaddress_[tNs[i]] = i;
        }
    }
    
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

                    toNeighbour << tNaddress_;
                }
            }
        }
    
        //- receiving
        for (int p = 0; p < Pstream::nProcs(); p++)
        {
            if(p != Pstream::myProcNo())
            {
                List<label> tNaddressProc;

                const int proc = p;
                {
                    IPstream fromNeighbour(Pstream::blocking, proc);
                    fromNeighbour  >> tNaddressProc;
                }
    
                forAll(tNaddressProc, i)
                {
                    if(tNaddressProc[i] != -1)
                    {
                        tNaddress_[i]=tNaddressProc[i];
                    }
                }
            }
        }
    }
    
    nLiquidMols_ = 0;
    
    forAll(tNaddress_, i)
    {
        if(tNaddress_[i] != -1)
        {
            nLiquidMols_++;
        }
    }
    
    Pout << "number of liquid mols = " << nLiquidMols_ << endl;
    
    
    
    acf_.setSize(nLiquidMols_);
    
    velocities_.clear();
    velocities_.setSize(nLiquidMols_);
    
    forAll(velocities_, i)
    {
        velocities_[i].setSize(nSteps_, 0.0);
        acf_[i].setSize(nSteps_, 0.0);
    }

    ACF_.setSize(nSteps_, 0.0);
    
    // set initial velocities 
        
    setVelocities();
    nS_++;
    
    Pout << "done " << endl;

}


void axialSelfDiffusion::setVelocities()
{
    IDLList<polyMolecule>::iterator mol(molCloud_.begin());
    
    for (mol = molCloud_.begin(); mol != molCloud_.end(); ++mol)
    {
        if(findIndex(molIds_, mol().id()) != -1)
        {
            label tN = mol().trackingNumber();
            label molIndex=tNaddress_[tN];
            
            if(molIndex != -1)
            {
                velocities_[molIndex][nS_] = mol().v() & unitVector_;
            }
        }
    }     
}

void axialSelfDiffusion::calculateField()
{
    if(nS_ > nSteps_)
    {
        Info << "axialSelfDiffusion averaging " << endl;
        
        // parallel processing
        
        if(Pstream::parRun())
        {
            forAll(velocities_, i)
            {
                forAll(velocities_[i], j)
                {
                    reduce(velocities_[i][j], sumOp<scalar>());
                }
            }
        } 
    
        // calculate molecule based ACF 
        
        label T = nSteps_;
        
        for (label i=0; i<nLiquidMols_; i++)
        {
            for (label k=0; k<T; k++)
            {
                scalar uSum = 0.0;
                
                for (label t=0; t<T-k; t++)
                {
                    uSum += velocities_[i][t]*velocities_[i][t+k];
                }
            
                acf_[i][k] += uSum/(T-1);
            }
        }
    

        nBatch_ += 1.0;        
        
        // the acf summed across molecule batches and all molecules
        
        ACF_ = 0.0; // reset
        
        for (label i=0; i<nLiquidMols_; i++)
        {
            for (label k=0; k<T; k++)
            {
                ACF_[k] += acf_[i][k]/(nBatch_*scalar(nLiquidMols_));
            }
        }
        
        
        // integrate ACF_ 
        
        
        // reset 
        
        nS_ = 0;
        
        forAll(velocities_, i)
        {
            forAll(velocities_[i], j)
            {
                velocities_[i][j] = 0.0;
            }
        }
    }
    else
    {
        setVelocities();
        nS_++;
    }    
}

// scalar axialSelfDiffusion::getIntegral()
// {
//     scalar timeIntegration = 0.0;
// 
//     const scalar& dt = time_.deltaT().value();
//     
//     if(((f.size() -1) % 2) == 0)// simpsons 1/3 rule
//     {
//         timeIntegration += f[0];
//         timeIntegration += f[f.size()-1];
// 
//         for (label i=1; i<f.size()-1; i++)
//         {
//             if((i % 2) == 0) // -even
//             {
//                 timeIntegration += 2.0*f[i];
//             }
//             else // odd
//             {
//                 timeIntegration += 4.0*f[i];
//             }
//         }
//         
//         timeIntegration *= dt/3.0;
// 
//     }
//     else // trapezoid rule
//     {
//         timeIntegration += f[0];
//         timeIntegration += f[f.size()-1];
// 
//         for (label i=1; i<f.size()-1; i++)
//         {
//             timeIntegration += 2.0*f[i];
//         }
// 
//         timeIntegration *= 0.5*dt;
//     }
//         
//     return timeIntegration;
// }

void axialSelfDiffusion::writeField()
{
    const Time& runTime = time_.time();

    if(runTime.outputTime())
    {
        if(Pstream::master())
        {
            scalarField timeFieldII (nSteps_, 0.0);
            scalarField acf (nSteps_, 0.0);
            
            const scalar& deltaT = time_.time().deltaT().value();
            
            forAll(timeFieldII, i)
            {
                timeFieldII[i]= deltaT*i;
                acf[i]=ACF_[i];
            }
            
            writeTimeData
            (
                timePath_,
                "ACF_"+fieldName_+".xy",
                timeFieldII,
                acf
            );
        }
    }
}

void axialSelfDiffusion::measureDuringForceComputation
(
    polyMolecule* molI,
    polyMolecule* molJ
){}

void axialSelfDiffusion::measureDuringForceComputationSite
(
    polyMolecule* molI,
    polyMolecule* molJ,
    label sI,
    label sJ
){}


const propertyField& axialSelfDiffusion::fields() const
{
    return fields_;
}



} // End namespace Foam

// ************************************************************************* //
