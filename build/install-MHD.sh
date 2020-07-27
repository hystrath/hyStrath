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

mkdir -p $sendingDir

# copy new files --------------------------------------------------------------
cp -r $currentDir/src $sendingDir/
cp -r $currentDir/applications $sendingDir/
cp -r $currentDir/run $sendingDir/


# compile new libraries -------------------------------------------------------
cd $sendingDir/src/mhdModels/
wclean libso
wmake -j$nProcs libso


# compile new executables ------------------------------------------------------
#---- solver ----
cd $sendingDir/applications/solvers/compressible/hy2MhdFoam
wclean
wmake -j$nProcs


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
Hybrid module CFD-MHD $WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
