#!/bin/bash

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$HOME/$WM_PROJECT/$userName-$WM_PROJECT_VERSION"

mkdir -p $sendingDir


# copy new files --------------------------------------------------------------
cp -r $currentDir/src $sendingDir/
cp -r $currentDir/applications $sendingDir/
cp -r $currentDir/run $sendingDir/

# compile new libraries -------------------------------------------------------
cd $sendingDir/src/thermophysicalModels/strath/
wclean all
./Allwmake

cd $sendingDir/src/TurbulenceModels/
wclean all
cd $sendingDir/src/TurbulenceModels/compressible
wmake libso
cd $sendingDir/src/TurbulenceModels/schemes
wmake libso

cd $sendingDir/src/thermophysicalModels/strath/
./AllwmakeBis

cd $sendingDir/src/hTCModels
wclean libso
wmake libso

cd $sendingDir/src/finiteVolume
wclean libso
wmake libso

cd $sendingDir/src/functionObjects/forces
wclean libso
cd $sendingDir/src/functionObjects
./Allwmake-hyStrath

cd $sendingDir/src/fvOptions
wclean libso
wmake libso


# compile new executables ------------------------------------------------------
#---- solvers ----
cd $sendingDir/applications/solvers/compressible/hy2Foam/
./Allwclean
./Allwmake

#---- utilities ----
cd $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
wclean
wmake

cd $sendingDir/applications/utilities/mesh/generation/blockMeshDG
wclean all
./Allwmake

cd $sendingDir/applications/utilities/postProcessing/wall/
./wcleanAll
./wmakeAll


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo -e "
cfdStrath_$WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
