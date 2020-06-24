#!/bin/bash

set -e

nProcs=1
if [ $# -ne 0 ]
  then nProcs=$1;
fi

tref_CFD=1020
tref_DSMC=1380
tref_PICDSMC=310
tref_MHD=160

display_bar() {
    local w=20 p=$1; shift
    local pct=$(( $p*$w/100 ))
    printf -v arrows "%*s" "$pct" ""; arrows=${arrows// />};
    printf "\r\e[K[%-*s] %3d%% %s" "$w" "$arrows" "$p" "$*"; 
}

progress_bar() {
    local slp=$1;
    for x in {1..99}
    do
        display_bar "$x" $2
        sleep $slp
    done ;
}

print_status() {
    local logfile="$1"
    local last_line=$( tail -n 2 $logfile )
    if [[ "$last_line" == *"successfully"* ]]; then
       echo -e "\nSUCCESS";
    else
       echo -e "\nFAIL: check $logfile";
    fi
}

install_CFD() {
    local t_CFD=$(($tref_CFD / $nProcs))
    local sleep_period=`bc -l <<< $t_CFD/100`
    progress_bar $sleep_period "installing CFD module" &
    progress_bar_pid=$!
    disown
    ./build/install-CFD.sh $nProcs > logInstall-CFD 2>&1
    display_bar "100" "installed CFD module"
    kill $progress_bar_pid >/dev/null 2>&1
    print_status logInstall-CFD
}

sync_CFD() {
    local t_CFD=$(($tref_CFD / $nProcs))
    local sleep_period=`bc -l <<< $t_CFD/100`
    progress_bar $sleep_period "syncing CFD module" &
    progress_bar_pid=$!
    disown
    ./build/resync-CFD.sh $nProcs > logSync-CFD 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "synced CFD module"
    print_status logSync-CFD
}

install_DSMC() {
    local t_DSMC=$(($tref_DSMC / $nProcs))
    local sleep_period=`bc -l <<< $t_DSMC/100`
    progress_bar $sleep_period "installing DSMC module" &
    progress_bar_pid=$!
    disown
    ./build/install-DSMC.sh $nProcs > logInstall-DSMC 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "installed DSMC module"
    print_status logInstall-DSMC
}

sync_DSMC() {
    local t_DSMC=$(($tref_DSMC / $nProcs))
    local sleep_period=`bc -l <<< $t_DSMC/100`
    progress_bar $sleep_period "syncing DSMC module" &
    progress_bar_pid=$!
    disown
    ./build/resync-DSMC.sh $nProcs > logSync-DSMC 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "synced DSMC module"
    print_status logSync-DSMC
}

install_PICDSMC() {
    local t_PICDSMC=$(($tref_PICDSMC / $nProcs))
    local sleep_period=`bc -l <<< $t_PICDSMC/100`
    progress_bar $sleep_period "installing hybrid PIC-DSMC module" &
    progress_bar_pid=$!
    disown
    ./build/install-hybridPICDSMC.sh $nProcs > logInstall-hybridPICDSMC 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "installed hybrid PIC-DSMC module"
    print_status logInstall-hybridPICDSMC
}

sync_PICDSMC() {
    local t_PICDSMC=$(($tref_PICDSMC / $nProcs))
    local sleep_period=`bc -l <<< $t_PICDSMC/100`
    progress_bar $sleep_period "syncing hybrid PIC-DSMC module" &
    progress_bar_pid=$!
    disown
    ./build/resync-hybridPICDSMC.sh $nProcs > logSync-hybridPICDSMC 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "synced hybrid PIC-DSMC module"
    print_status logSync-hybridPICDSMC
}

install_MHD() {
    local t_MHD=$(($tref_MHD / $nProcs))
    local sleep_period=`bc -l <<< $t_MHD/100`
    progress_bar $sleep_period "installing CFD-MHD module" &
    progress_bar_pid=$!
    disown
    ./build/install-MHD.sh $nProcs > logInstall-MHD 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "installed CFD-MHD module"
    print_status logInstall-MHD
}

sync_MHD() {
    local t_MHD=$(($tref_MHD / $nProcs))
    local sleep_period=`bc -l <<< $t_MHD/100`
    progress_bar $sleep_period "syncing CFD-MHD module" &
    progress_bar_pid=$!
    disown
    ./build/resync-MHD.sh $nProcs > logSync-MHD 2>&1
    kill $progress_bar_pid >/dev/null 2>&1
    display_bar "100" "synced CFD-MHD module"
    print_status logSync-MHD
}

echo "-----  INSTALLATION  -----"
echo "1 - CFD module"
echo "2 - DSMC module"
echo "3 - Hybrid PIC-DSMC module"
echo "4 - CFD-MHD module"
echo -e "5 - All modules\n"

echo "-----  SYNCHRONISATION  -----"
echo "11 - CFD module"
echo "12 - DSMC module"
echo "13 - Hybrid PIC-DSMC module"
echo "14 - CFD-MHD module"
echo -e "15 - All modules\n"

read -r -p "Enter choice: " input
 
case $input in
    1)
        install_CFD
        ;;
    2)
        install_DSMC
        ;;
    3)
        install_PICDSMC
        ;;
    4)
        install_MHD
        ;;
    5)
        install_CFD
        install_DSMC
        install_PICDSMC
        install_MHD
        ;;
    11)
        sync_CFD
        ;;
    12)
        sync_DSMC
        ;;
    13)
        sync_PICDSMC
        ;;
    14)
        sync_MHD
        ;;
    15)
        sync_CFD
        sync_DSMC
        sync_PICDSMC
        sync_MHD
        ;;          
    *)
 echo "Invalid input"
 exit 1
 ;;
esac
