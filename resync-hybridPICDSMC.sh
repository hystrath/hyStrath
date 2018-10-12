#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$WM_PROJECT_USER_DIR"

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi

mkdir -p $sendingDir

# synchronize folders and compile new libraries -------------------------------
rsync -rtvuc $currentDir/src/lagrangian/regionsAndInterfaces/ $sendingDir/src/lagrangian/regionsAndInterfaces/
cd $sendingDir/src/lagrangian/regionsAndInterfaces/
wmake -j$nProcs libso

rsync -rtvuc $currentDir/src/lagrangian/pic/ $sendingDir/src/lagrangian/pic/
cd $sendingDir/src/lagrangian/pic/
wmake -j$nProcs libso


# compile new executables ------------------------------------------------------
#---- solver ----
rsync -rtvuc $currentDir/applications/solvers/hybridMethods/pdFoam/ $sendingDir/applications/solvers/hybridMethods/pdFoam/
cd $sendingDir/applications/solvers/hybridMethods/pdFoam
wmake -j$nProcs

#---- utilities ----
rsync -rtvuc $currentDir/applications/utilities/preProcessing/hybrid/pdInitialise/ $sendingDir/applications/utilities/preProcessing/hybrid/pdInitialise/
cd $sendingDir/applications/utilities/preProcessing/hybrid/pdInitialise/
wmake -j$nProcs all


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
Hybrid module PIC-DSMC $WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
