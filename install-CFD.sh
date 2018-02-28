#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$HOME/$WM_PROJECT/$userName-$WM_PROJECT_VERSION"

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
cd $sendingDir/src/thermophysicalModels/strath/
wclean all
./Allwmake -j$nProcs

cd $sendingDir/src/TurbulenceModels/
wclean all
cd $sendingDir/src/TurbulenceModels/compressible
wmake -j$nProcs libso
cd $sendingDir/src/TurbulenceModels/schemes
wmake -j$nProcs libso

cd $sendingDir/src/thermophysicalModels/strath/
./AllwmakeBis -j$nProcs

cd $sendingDir/src/hTCModels
wclean libso
wmake -j$nProcs libso

cd $sendingDir/src/finiteVolume
wclean libso
wmake -j$nProcs libso

cd $sendingDir/src/functionObjects/forces
wclean libso
cd $sendingDir/src/functionObjects
./Allwmake-hyStrath -j$nProcs

cd $sendingDir/src/fvOptions
wclean libso
wmake -j$nProcs libso


# compile new executables ------------------------------------------------------
#---- solvers ----
cd $sendingDir/applications/solvers/compressible/hy2Foam/
./Allwclean
./Allwmake -j$nProcs

#---- utilities ----
cd $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
wclean
wmake -j$nProcs

cd $sendingDir/applications/utilities/mesh/generation/blockMeshDG
wclean all
./Allwmake -j$nProcs

cd $sendingDir/applications/utilities/postProcessing/wall/
./wcleanAll
./wmakeAll -j$nProcs


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
cfdStrath_$WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
