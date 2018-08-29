#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

set -e

userName=`whoami`

currentDir=`pwd`
sendingDir="$FOAM_INST_DIR/$userName-$WM_PROJECT_VERSION"

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi

mkdir -p $sendingDir


# copy new files --------------------------------------------------------------
foldersSrc="lagrangian parallel functionObjects/field-dsmcStrath functionObjects/lagrangian"
filesInFolderSrc="functionObjects"
foldersApp="solvers/discreteMethods utilities/preProcessing/mapFields utilities/preProcessing/mapFieldsPar utilities/preProcessing/dsmc utilities/parallelProcessing utilities/postProcessing/miscellaneous/foamListTimes utilities/mesh/generation/makeAxialMesh utilities/mesh/generation/blockMeshDG"

for folder in $foldersSrc
do
  mkdir -p $sendingDir/src/$folder
  cp -r $currentDir/src/$folder $sendingDir/src/`dirname $folder`
done

for filesInFolder in $filesInFolderSrc
do
  find $currentDir/src/$filesInFolder/ -maxdepth 1 -type f | xargs cp -t $sendingDir/src/$filesInFolder
done

for folder in $foldersApp
do
  mkdir -p $sendingDir/applications/$folder
  cp -r $currentDir/applications/$folder $sendingDir/applications/`dirname $folder`
done

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

cd $sendingDir/src/functionObjects/field-dsmcStrath
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
DSMC module $WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
