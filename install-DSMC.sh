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
cd $sendingDir/src/lagrangian/
wclean all
./Allwmake -j$nProcs

cd $sendingDir/src/parallel/decompose/decompose/
wclean libso
wmake -j$nProcs libso

cd $sendingDir/src/parallel/reconstruct/
wclean all
./Allwmake -j$nProcs

cd $sendingDir/src/functionObjects/field
wclean libso
cd $sendingDir/src/functionObjects/lagrangian
wclean libso
cd $sendingDir/src/functionObjects/
./Allwmake-dsmcStrath -j$nProcs


# compile new executables ------------------------------------------------------
#---- solvers ----
cd $sendingDir/applications/solvers/discreteMethods/dsmc/
./wcleanAll
./wmakeAll -j$nProcs

cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdEquilibrationFoam/
wclean
wmake -j$nProcs

cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdFoam/
wclean
wmake -j$nProcs

#---- utilities ----
cd $sendingDir/applications/utilities/preProcessing/mapFields
wclean
wmake -j$nProcs

cd $sendingDir/applications/utilities/preProcessing/mapFieldsPar
wclean
wmake -j$nProcs

cd $sendingDir/applications/utilities/preProcessing/dsmc/
./wcleanAll
./wmakeAll -j$nProcs

cd $sendingDir/applications/utilities/parallelProcessing/
wclean all
wmake  -j$nProcs all

cd $sendingDir/applications/utilities/postProcessing/miscellaneous/foamListTimes/
wclean
wmake -j$nProcs

cd $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
wclean
wmake -j$nProcs

cd $sendingDir/applications/utilities/mesh/generation/blockMeshDG
wclean all
./Allwmake -j$nProcs


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
dsmcStrath_$WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
