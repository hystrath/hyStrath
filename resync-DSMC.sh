#!/bin/bash

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$HOME/$WM_PROJECT/$userName-$WM_PROJECT_VERSION"

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi


# synchronize folders and compile new libraries -------------------------------
rsync -rtvuc $currentDir/src/lagrangian/ $sendingDir/src/lagrangian/
cd $sendingDir/src/lagrangian/
./Allwmake -j$nProcs

rsync -rtvuc $currentDir/src/parallel/decompose/decompose/ $sendingDir/src/parallel/decompose/decompose/
cd $sendingDir/src/parallel/decompose/decompose/
wmake -j$nProcs libso

rsync -rtvuc $currentDir/src/parallel/reconstruct/ $sendingDir/src/parallel/reconstruct/
cd $sendingDir/src/parallel/reconstruct/
./Allwmake -j$nProcs

rsync -rtvuc $currentDir/src/functionObjects/field $sendingDir/src/functionObjects/field
rsync -rtvuc $currentDir/src/functionObjects/lagrangian $sendingDir/src/functionObjects/lagrangian
cd $sendingDir/src/functionObjects/
./Allwmake-dsmcStrath -j$nProcs


# compile new executables ------------------------------------------------------
#---- solvers ----
rsync -rtvuc $currentDir/applications/solvers/discreteMethods/dsmc/ $sendingDir/applications/solvers/discreteMethods/dsmc/
cd $sendingDir/applications/solvers/discreteMethods/dsmc/
./wmakeAll -j$nProcs

rsync -rtvuc $currentDir/applications/solvers/discreteMethods/molecularDynamics/mdEquilibrationFoam/ $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdEquilibrationFoam/
cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdEquilibrationFoam/
wmake -j$nProcs

rsync -rtvuc $currentDir/applications/solvers/discreteMethods/molecularDynamics/mdFoam/ $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdFoam/
cd $sendingDir/applications/solvers/discreteMethods/molecularDynamics/mdFoam/
wmake -j$nProcs

#---- utilities ----
rsync -rtvuc $currentDir/applications/utilities/preProcessing/mapFields $sendingDir/applications/utilities/preProcessing/mapFields
cd $sendingDir/applications/utilities/preProcessing/mapFields
wmake -j$nProcs

rsync -rtvuc $currentDir/applications/utilities/preProcessing/mapFieldsPar $sendingDir/applications/utilities/preProcessing/mapFieldsPar
cd $sendingDir/applications/utilities/preProcessing/mapFieldsPar
wmake -j$nProcs

rsync -rtvuc $currentDir/applications/utilities/preProcessing/dsmc/ $sendingDir/applications/utilities/preProcessing/dsmc/
cd $sendingDir/applications/utilities/preProcessing/dsmc/
./wmakeAll -j$nProcs

rsync -rtvuc $currentDir/applications/utilities/parallelProcessing/ $sendingDir/applications/utilities/parallelProcessing/
cd $sendingDir/applications/utilities/parallelProcessing/
wmake  -j$nProcs all

rsync -rtvuc $currentDir/applications/utilities/postProcessing/miscellaneous/foamListTimes/ $sendingDir/applications/utilities/postProcessing/miscellaneous/foamListTimes/
cd $sendingDir/applications/utilities/postProcessing/miscellaneous/foamListTimes/
wmake -j$nProcs

rsync -rtvuc $currentDir/applications/utilities/mesh/generation/makeAxialMesh $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
cd $sendingDir/applications/utilities/mesh/generation/makeAxialMesh
wmake -j$nProcs

rsync -rtvuc $currentDir/applications/utilities/mesh/generation/blockMeshDG $sendingDir/applications/utilities/mesh/generation/blockMeshDG
cd $sendingDir/applications/utilities/mesh/generation/blockMeshDG
./Allwmake -j$nProcs



# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
dsmcStrath_$WM_PROJECT_VERSION updated successfully. Hope you'll enjoy it, $userName :)
"
