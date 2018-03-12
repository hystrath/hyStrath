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
foldersSrc="thermophysicalModels TurbulenceModels hTCModels finiteVolume fvOptions functionObjects/forces"
filesInFolderSrc="functionObjects"
foldersApp="solvers/compressible/hy2Foam utilities/mesh/generation/makeAxialMesh utilities/mesh/generation/blockMeshDG"

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
./Allwmake-cfdStrath -j$nProcs

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


# re-set to the initial directory ---------------------------------------------
cd $currentDir

echo "
cfdStrath_$WM_PROJECT_VERSION compiled successfully. Hope you'll enjoy it, $userName :)
"
