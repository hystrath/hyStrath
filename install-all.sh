#!/bin/bash

set -e

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi

./install-CFD.sh $nProcs
./install-DSMC.sh $nProcs
./install-hybridPICDSMC.sh $nProcs
./install-MHD.sh $nProcs
