/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Class
    dsmcDynamicLoadBalancing

Description

\*----------------------------------------------------------------------------*/

#include "dsmcDynamicLoadBalancing.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



//- Constructor
dsmcDynamicLoadBalancing::dsmcDynamicLoadBalancing
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    dsmcLoadBalanceDict_
    (
        IOobject
        (
            "loadBalanceDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    /*controlDict_
    (
        IOobject
        (
            "controlDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),*/
    performBalance_(false),
    enableBalancing_(Switch(dsmcLoadBalanceDict_.lookup("enableBalancing"))),
    originalEndTime_(time_.time().endTime().value()),
    maxImbalance_(readScalar(dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")))
    //nProcs_(readLabel(dsmcLoadBalanceDict_.lookup("numberOfSubdomains")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::update()
{
    if(time_.time().outputTime())
    {
        //- Checking for modifications in the IOdictionary
        //  this allows for run-time tuning of any parameters.  
        // DLETED VINCENT: not useful
        
        /*IOdictionary newDict
        (
            IOobject
            (
                "loadBalanceDict",
                time_.system(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );*/

        updateProperties(/*newDict*/);
        
        //- Load Balancing
        if (Pstream::parRun())
        {
            const scalar& allowableImbalance = maxImbalance_;
                
            // First determine current level of imbalance - do this for all
            // parallel runs, even if balancing is disabled
            scalar nGlobalParticles = cloud_.size();
            Foam::reduce(nGlobalParticles, sumOp<scalar>());
            
            scalar idealNParticles = scalar(nGlobalParticles)/scalar(Pstream::nProcs());
            
            scalar nParticles = cloud_.size();
            scalar localImbalance = mag(nParticles - idealNParticles);
            Foam::reduce(localImbalance, maxOp<scalar>());
            scalar maxImbalance = localImbalance/idealNParticles;
            
            Info << "    Maximum imbalance = " << 100*maxImbalance << "%" << endl;
            
            if( maxImbalance > allowableImbalance && enableBalancing_)
            {   
                performBalance_ = true;
                
                originalEndTime_ = time_.time().endTime().value();
                
                scalar currentTime = time_.time().value();
                
                time_.setEndTime(currentTime);
            }
        }
    }
}
    
    
void dsmcDynamicLoadBalancing::copyPolyMeshToLatestTimeFolder() const
{
    const fileName constantInProcessor0 = "processor0/constant";
    
    if (not isDir(constantInProcessor0))
    {
        for (label i=0; i<Pstream::nProcs(); i++)
        {
            const word findStartTime = 
                "starting=`foamListTimes -processor -withZero -startTime`; ";
                
            if (findStartTime != "0")
            {    
                const word findLatestTime = 
                    "latest=`foamListTimes -processor -withZero -latestTime`; ";
                const word findTimes = findStartTime + findLatestTime;
                const word processorName = "processor" + name(i) + "/";
                const word copyPolyMesh = findTimes + "cp -r " 
                    + processorName + "$starting" + "/polyMesh " 
                    + processorName + "$latest/";
                
                Foam::system(copyPolyMesh);
            }
        }
    }
    /*else
    {
        const word findLatestTime = 
                    "latest=`foamListTimes -processor -withZero -latestTime`; ";
                    
        const fileName polyMeshInTimeFolder = "processor0/" + findLatestTime + "/polyMesh";
        
        if (isDir(polyMeshInTimeFolder))
        {
            const fileName constantInProcessor0 = "processor0/constant";
        }
    }*/
}


void dsmcDynamicLoadBalancing::perform(const int noRefinement)
{    
    if(performBalance_)
    {      
        if (Pstream::master())
        {   
            if (noRefinement == 0)
            {
                copyPolyMeshToLatestTimeFolder();
            }
            
            //string redistributeCommand("mpirun -np " + word(nProcs_) + " redistributeParDSMCLoadBalance -parallel");
            //system("redistributeCommand");

            /*const char reconstructParMeshCommand[] =  
                "reconstructParMesh -latestTime";
                
            const word reconstructParCommand = 
                "reconstructPar -latestTime -parallel"; 
                
            const word decomposeDSMCLoadBalanceParCommand = 
                "decomposeDSMCLoadBalancePar -force -latestTime -copyUniform -parallel";       
            
            const word decomposeDSMCCommand = word("var=`foamListTimes -noZero`; ")
                + word("if [ -z ${var+x} ]; then echo 'error'; else ") 
                + decomposeDSMCLoadBalanceParCommand + word("; fi");*/
                
            // Open MPI does not support recursive calls of mpirun
            // MPI_Comm_spawn must be used instead
            //int threading_level_required = MPI_THREAD_MULTIPLE;
            //int threading_level_provided;
            //MPI_Init_thread(0, 0, threading_level_required, &threading_level_provided);
            /*MPI_Comm comm_to_workers;
            char worker_program[100];
            strcpy(worker_program, "hhh");*/
            
            /*MPI_Comm_spawn(worker_program, MPI_ARGV_NULL, Pstream::nProcs(),
                MPI_INFO_NULL, 0, MPI_COMM_SELF, &comm_to_workers, MPI_ERRCODES_IGNORE);*/
                
            //system(reconstructParMeshCommand);
            //system(reconstructParCommand);
            //system("rm -r processor*");
            //system(decomposeDSMCCommand);
            
            Foam::system("reconstructParMesh -latestTime");
            Foam::system("reconstructPar -latestTime");
            Foam::system("rm -r processor*");
            
            Foam::system("var=`foamListTimes -noZero`; if [ -z ${var+x} ]; then echo 'error'; else decomposeDSMCLoadBalancePar -force -latestTime -copyUniform; fi");
            
            // Backup folders must be stored in resultFolders and moved back
            // as the simulation finishes
            //system("foamListTimes -rm"); // no backup alternative but risky
            mkDir("resultFolders");
            Foam::system("var2=`foamListTimes`; mv $var2 resultFolders");
            
            performBalance_ = false;
        }
            
        time_.setEndTime(originalEndTime_);
    }
}

void dsmcDynamicLoadBalancing::updateProperties
(
    /*const IOdictionary& newDict*/
)
{
    enableBalancing_ = Switch(dsmcLoadBalanceDict_.lookup("enableBalancing"));
    maxImbalance_ = readScalar(dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance"));
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
