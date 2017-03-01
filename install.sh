#!/bin/bash

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$HOME/OpenFOAM/$userName-$WM_PROJECT_VERSION"

mkdir -p $sendingDir

# copy new files --------------------------------------------------------------
cp -r $currentDir/src $sendingDir/
cp -r $currentDir/applications $sendingDir/
cp -r $currentDir/run $sendingDir/

# compile new libraries -------------------------------------------------------
cd $sendingDir/src/thermophysicalModels/strath/
./Allwmake

cd $sendingDir/src/turbulenceModels/strath/
./Allwmake

cd $sendingDir/src/hTCModels
wclean libso
wmake libso

cd $sendingDir/src/finiteVolume
wclean libso
wmake libso

cd $sendingDir/src/postProcessing/functionObjects/forces
wclean libso
wmake libso


# compile new executables ------------------------------------------------------
#---- solvers ----
cd $sendingDir/applications/solvers/compressible/hy2Foam/
./Allwmake

#---- utilities ----
cd $sendingDir/applications/utilities/postProcessing/wall/pressureCoeff
wclean
wmake

cd $sendingDir/applications/utilities/postProcessing/wall/frictionCoeff2T
wclean
wmake

cd $sendingDir/applications/utilities/postProcessing/wall/yPlus2T
wclean
wmake

cd $sendingDir/applications/utilities/postProcessing/wall/StantonNumber
wclean
wmake




# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo -e "
cfdStrath has now been compiled. Hope you'll enjoy it $userName :)
"
