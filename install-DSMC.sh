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
cd $sendingDir/src/lagrangian/
wclean all
./Allwmake

cd $sendingDir/src/parallel/decompose/decompose/
wclean libso
wmake libso

cd $sendingDir/src/parallel/reconstruct/
wclean all
./Allwmake

cd $sendingDir/src/functionObjects/field
wclean libso
cd $sendingDir/src/functionObjects/lagrangian
wclean libso
cd $sendingDir/src/functionObjects/
./Allwmake-dsmcStrath


# compile new executables ------------------------------------------------------
#---- solvers ----
cd $sendingDir/applications/solvers/discreteMethods/dsmc/
./wcleanAll
./wmakeAll

cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdEquilibrationFoam/
wclean
wmake

cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdFoam/
wclean
wmake

#---- utilities ----
cd $sendingDir/applications/utilities/preProcessing/mapFields
wclean
wmake

cd $sendingDir/applications/utilities/preProcessing/mapFieldsPar
wclean
wmake

cd $sendingDir/applications/utilities/preProcessing/dsmc/
./wcleanAll
./wmakeAll

cd $sendingDir/applications/utilities/parallelProcessing/
wclean all
wmake all

cd $sendingDir/applications/utilities/postProcessing/miscellaneous/foamListTimes/
wclean
wmake

cd $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
wclean
wmake

cd $sendingDir/applications/utilities/mesh/generation/blockMeshDG
wclean all
./Allwmake


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo -e "
dsmcStrath_$WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
