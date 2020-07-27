#!/bin/bash
cd ${0%/*}/.. || exit 1    # Run from hyStrath directory

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$WM_PROJECT_USER_DIR"

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi


# synchronize folders and compile new libraries -------------------------------
rsync -rtvuc $currentDir/src/mhdModels/ $sendingDir/src/mhdModels/
cd $sendingDir/src/mhdModels/
wmake -j$nProcs libso


# compile new executables ------------------------------------------------------
#---- solvers ----
rsync -rtvuc $currentDir/applications/solvers/compressible/hy2MhdFoam/ $sendingDir/applications/solvers/compressible/hy2MhdFoam/
cd $sendingDir/applications/solvers/compressible/hy2MhdFoam/
wmake -j$nProcs


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
Hybrid module CFD-MHD $WM_PROJECT_VERSION updated successfully. Hope you'll enjoy it, $userName :)
"
