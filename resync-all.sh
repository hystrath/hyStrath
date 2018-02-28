#!/bin/bash

set -e

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi

./resync-CFD.sh $nProcs
./resync-DSMC.sh $nProcs
