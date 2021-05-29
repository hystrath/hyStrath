/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
    balanceUntilTime_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<scalar>
        (
            "balanceUntilTime",
            VGREAT
        )
    ),
    originalEndTime_(time_.time().endTime().value()),
    maxImbalance_(readScalar
    (
        dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")
    )),
    limitTimeDirBackups_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "limitTimeDirBackups",
            -1
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::update()
{
    if (time_.time().outputTime())
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

            Info<< "    Maximum imbalance = " << 100*maxImbalance << "%" << nl
                << endl;

            // explanation of modes:
            // 1. enableBalancing = true:
            //   1. if time <= balanceUntilTime -> load balance
            //   2. if time > balanceUntilTime -> do not load balance
            // 2. enableBalancing = false -> do not load balance
            if
            (
                   enableBalancing_
                && time_.time().value() <= balanceUntilTime_
                && maxImbalance > allowableImbalance
            )
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
                "starting=$(foamListTimes -processor -withZero -startTime);";

            if (findStartTime != "0")
            {
                const word findLatestTime =
                    "latest=$(foamListTimes -processor -withZero -latestTime);";
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
    if (performBalance_)
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
                word("timeDirs=$(foamListTimes -noZero);")
                // check if there are any time dirs, if not this indicates a
                // fatal error
                + word("if [ -z ${timeDirs+x} ];")
                + word("then echo \"error\";")
                // decompose the latest time dir
                + word("else decomposeDSMCLoadBalancePar ")
                + word("-force -latestTime -copyUniform;")
                + word("fi");
            Foam::system(decomposePar);

            // backup time dirs must be stored in resultFolders to prevent them
            // from being cleared. They can be moved back when the simulation
            // has finished.
            mkDir("resultFolders");

            // respect the specified limit of max. concurrent time dir backups
            if (limitTimeDirBackups_ >= 0)
            {
                // impose the limit by moving back the time dirs that are
                // currently already backuped up and then only keeping the
                // most recent time dirs <= the imposed limit.
                const word backupTimeDirsWithLimit =
                    // move the time dirs currently backed up to the case dir
                    // so the foamListTimes utility can be used
                    // make sure directory is not empty to prevent mv from
                    // printing a warning
                    word("if [ \"$(ls -A resultFolders)\" ];")
                    + word("then mv resultFolders/* .; fi;")
                    // total number of time directories
                    + word("timeDirs=$(foamListTimes);")
                    + word("nTimeDirs=$(echo $timeDirs | tr -cd ' ' | wc -c);")
                    + word("nTimeDirs=$((nTimeDirs+1));")
                    // convert OpenFOAM label to shell variable for limit
                    + word("limitNTimeDirs=")
                    + name(limitTimeDirBackups_)
                    + word(";")
                    // if the total number of time directories is larger than
                    // limit remove the oldest time directories first
                    + word("diffNTimeDirs=$((nTimeDirs-limitNTimeDirs));")
                    + word("if [ $diffNTimeDirs -gt 0 ];")
                    + word("then timeDirs=$(echo $timeDirs | ")
                    + word("cut -d ' ' -f$((diffNTimeDirs+1))-$nTimeDirs);")
                    + word("fi;")
                    // move the time dirs that were chosen to be eligible for
                    // backup to the backup dir
                    + word("mv $timeDirs resultFolders/;")
                    // clear all other time dirs
                    + word("foamListTimes -rm");
                Foam::system(backupTimeDirsWithLimit);
            }
            else
            {
                // keep all time dirs in the backup dir
                Foam::system
                (
                    "timeDirs=$(foamListTimes); mv $timeDirs resultFolders/"
                );
            }

            performBalance_ = false;
        }

        time_.setEndTime(originalEndTime_);
    }
}

void dsmcDynamicLoadBalancing::updateProperties()
{
    enableBalancing_ = Switch(dsmcLoadBalanceDict_.lookup("enableBalancing"));
    // if balancing is active this additional option allows to specify a time
    // after which balancing is deactivated (this is useful in conjunction with
    // resetAtOutput / resetAtOutputUntilTime and averaging across solver
    // restarts)
    balanceUntilTime_ = dsmcLoadBalanceDict_.lookupOrDefault<scalar>
    (
        "balanceUntilTime",
        VGREAT
    );
    maxImbalance_ = readScalar
    (
        dsmcLoadBalanceDict_.lookup("maximumAllowableImbalance")
    );
    limitTimeDirBackups_ = dsmcLoadBalanceDict_.lookupOrDefault<label>
    (
        "limitTimeDirBackups",
        -1
    );
}


}  // End namespace Foam

// ************************************************************************* //
