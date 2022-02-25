#!/bin/bash

nProcs=1
if [ $# -ne 0 ]
then
    nProcs=$1
fi

tref_CFD=1200
tref_DSMC=1380
tref_PICDSMC=310

speedup=1.0
if [[ $nProcs -eq 2 || $nProcs -eq 3 ]]
then
    speedup=1.5
elif [ $nProcs -ge 4 ]; then
    speedup=2.0
fi

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
    done
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
    local sleep_period=`bc -l <<< $tref_CFD/$speedup/100`
    progress_bar $sleep_period "installing CFD module" &
    ./build/install-CFD.sh $nProcs > logInstall-CFD 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "installed CFD module"
    print_status logInstall-CFD
}

sync_CFD() {
    local sleep_period=`bc -l <<< $tref_CFD/$speedup/100`
    progress_bar $sleep_period "syncing CFD module" &
    ./build/resync-CFD.sh $nProcs > logSync-CFD 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "synced CFD module"
    print_status logSync-CFD
}

install_DSMC() {
    local sleep_period=`bc -l <<< $tref_DSMC/$speedup/100`
    progress_bar $sleep_period "installing DSMC module" &
    ./build/install-DSMC.sh $nProcs > logInstall-DSMC 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "installed DSMC module"
    print_status logInstall-DSMC
}

sync_DSMC() {
    local sleep_period=`bc -l <<< $tref_DSMC/$speedup/100`
    progress_bar $sleep_period "syncing DSMC module" &
    ./build/resync-DSMC.sh $nProcs > logSync-DSMC 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "synced DSMC module"
    print_status logSync-DSMC
}

install_PICDSMC() {
    local sleep_period=`bc -l <<< $tref_PICDSMC/$speedup/100`
    progress_bar $sleep_period "installing hybrid PIC-DSMC module" &
    ./build/install-hybridPICDSMC.sh $nProcs > logInstall-hybridPICDSMC 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "installed hybrid PIC-DSMC module"
    print_status logInstall-hybridPICDSMC
}

sync_PICDSMC() {
    local sleep_period=`bc -l <<< $tref_PICDSMC/$speedup/100`
    progress_bar $sleep_period "syncing hybrid PIC-DSMC module" &
    ./build/resync-hybridPICDSMC.sh $nProcs > logSync-hybridPICDSMC 2>&1 &
    wait %2
    if [ ! -z `jobs -p` ]
    then
      progress_bar_pid=$!
      disown
      pkill -P $$  > /dev/null 2>&1
    fi
    display_bar "100" "synced hybrid PIC-DSMC module"
    print_status logSync-hybridPICDSMC
}

echo "-----  INSTALLATION  -----"
echo "1 - CFD module"
echo "2 - DSMC module"
echo "3 - Hybrid PIC-DSMC module"
echo -e "4 - All modules\n"

echo "-----  SYNCHRONISATION  -----"
echo "11 - CFD module"
echo "12 - DSMC module"
echo "13 - Hybrid PIC-DSMC module"
echo -e "14 - All modules\n"

echo -e "Enter choice: \c"
read input
 
case $input in
    1)
        install_CFD
        ;;
    2)
        install_DSMC
        ;;
    3)
        echo -e "Confirm that the DSMC module is already installed [Y/n]: \c"
        read input_conf
        case $input_conf in
            [yY][eE][sS]|[yY])
                install_PICDSMC
                ;;
            [nN][oO]|[nN])
                echo "Please install the DSMC module first"
                ;;
            *)
                echo "Invalid input"
                exit 1
                ;;
        esac
        ;;
    4)
        install_CFD
        wait -n
        install_DSMC
        wait -n
        install_PICDSMC
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
    15)
        sync_CFD
        wait -n
        sync_DSMC
        wait -n
        sync_PICDSMC
        ;;             
    *)
 echo "Invalid input"
 exit 1
 ;;
esac
