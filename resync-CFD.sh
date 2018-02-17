#!/bin/bash

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$HOME/$WM_PROJECT/$userName-$WM_PROJECT_VERSION"

mkdir -p $sendingDir

# synchronize folders and compile new libraries -------------------------------
rsync -rtvuc $currentDir/src/thermophysicalModels/strath/ $sendingDir/src/thermophysicalModels/strath/
cd $sendingDir/src/thermophysicalModels/strath/
./Allwmake

rsync -rtvuc $currentDir/src/TurbulenceModels/compressible/ $sendingDir/src/TurbulenceModels/compressible/
cd $sendingDir/src/TurbulenceModels/compressible
wmake libso
rsync -rtvuc $currentDir/src/TurbulenceModels/schemes/ $sendingDir/src/TurbulenceModels/schemes/
cd $sendingDir/src/TurbulenceModels/schemes
wmake libso

cd $sendingDir/src/thermophysicalModels/strath/
./AllwmakeBis

rsync -rtvuc $currentDir/src/hTCModels/ $sendingDir/src/hTCModels/
cd $sendingDir/src/hTCModels
wmake libso

rsync -rtvuc $currentDir/src/finiteVolume/ $sendingDir/src/finiteVolume/
cd $sendingDir/src/finiteVolume
wmake libso

rsync -rtvuc $currentDir/src/functionObjects/ $sendingDir/src/functionObjects/
cd $sendingDir/src/functionObjects
./Allwmake-hyStrath

rsync -rtvuc $currentDir/src/fvOptions/ $sendingDir/src/fvOptions/
cd $sendingDir/src/fvOptions
wmake libso


# compile new executables ------------------------------------------------------
#---- solvers ----
rsync -rtvuc $currentDir/applications/solvers/compressible/hy2Foam/ $sendingDir/applications/solvers/compressible/hy2Foam/
cd $sendingDir/applications/solvers/compressible/hy2Foam/
./Allwmake

#---- utilities ----
rsync -rtvuc $currentDir/applications/utilities/generation/makeAxialMesh $sendingDir/applications/utilities/generation/makeAxialMesh
cd $sendingDir/applications/utilities/generation/makeAxialMesh
wmake

rsync -rtvuc $currentDir/applications/utilities/generation/blockMeshDG $sendingDir/applications/utilities/generation/blockMeshDG
cd $sendingDir/applications/utilities/generation/blockMeshDG
./Allwmake

rsync -rtvuc $currentDir/applications/utilities/postProcessing/wall $sendingDir/applications/utilities/postProcessing/wall
cd $sendingDir/applications/utilities/postProcessing/wall
./wmakeAll


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo -e "
cfdStrath_$WM_PROJECT_VERSION updated successfully. Hope you'll enjoy it, $userName :)
"
