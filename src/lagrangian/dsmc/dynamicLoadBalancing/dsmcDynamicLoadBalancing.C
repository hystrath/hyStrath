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
    performBalance_(false),
    enableBalancing_(Switch(dsmcLoadBalanceDict_.lookup("enableBalancing"))),
    originalEndTime_(time_.time().endTime().value()),
    maxImbalance_(readScalar
    (
        dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")
    ))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::update()
{
    if(time_.time().outputTime())
    {
        updateProperties();

        //- Load Balancing
        if (Pstream::parRun())
        {
            const scalar& allowableImbalance = maxImbalance_;

            // First determine current level of imbalance - do this for all
            // parallel runs, even if balancing is disabled
            scalar nGlobalParticles = cloud_.size();
            Foam::reduce(nGlobalParticles, sumOp<scalar>());

            scalar idealNParticles =
                scalar(nGlobalParticles)/scalar(Pstream::nProcs());

            scalar nParticles = cloud_.size();
            scalar localImbalance = mag(nParticles - idealNParticles);
            Foam::reduce(localImbalance, maxOp<scalar>());
            scalar maxImbalance = localImbalance/idealNParticles;

            Info<< "    Maximum imbalance = " << 100*maxImbalance << "%"
                << endl;

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
                "starting=`foamListTimes -processor -withZero -startTime`;";

            if (findStartTime != "0")
            {
                const word findLatestTime =
                    "latest=`foamListTimes -processor -withZero -latestTime`;";
                const word findTimes = findStartTime + findLatestTime;
                const word processorName = "processor" + name(i) + "/";
                const word copyPolyMesh = findTimes + "cp -r "
                    + processorName + "$starting" + "/polyMesh "
                    + processorName + "$latest/";

                Foam::system(copyPolyMesh);
            }
        }
    }
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
            Foam::system("reconstructParMesh -latestTime");
            Foam::system("reconstructPar -latestTime");
            Foam::system("rm -r processor*");

            const word decomposePar =
                word("timeDirs=`foamListTimes -noZero`;")
                // check if there are any time dirs, if not this indicates a
                // fatal error
                + word("if [ -z ${timeDirs+x} ];")
                + word("then echo \"error\";")
                // decompose the latest time dir
                + word("else decomposeDSMCLoadBalancePar ")
                + word("-force -latestTime -copyUniform;")
                + word("fi");
            Foam::system(decomposePar);

            // Backup folders must be stored in resultFolders and moved back
            // when the simulation finishes
            mkDir("resultFolders");
            Foam::system
            (
                "timeDirs=`foamListTimes`; mv $timeDirs resultFolders/"
            );

            performBalance_ = false;
        }

        time_.setEndTime(originalEndTime_);
    }
}

void dsmcDynamicLoadBalancing::updateProperties()
{
    enableBalancing_ = Switch(dsmcLoadBalanceDict_.lookup("enableBalancing"));
    maxImbalance_ = readScalar
    (
        dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")
    );
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

}  // End namespace Foam

// ************************************************************************* //
