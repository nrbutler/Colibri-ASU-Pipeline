#!/bin/bash

[ -f /home/colibri/MANUAL_MODE ] && { echo "file /home/colibri/MANUAL_MODE present, aborting" ; exit 1 ; }

source /usr/local/var/colibri/bin/redux_funcs_colibri.sh
test=

# just in case
rm $REDUX_LOCKFILE /home/colibri/monitoring.txt 2>/dev/null

# perhaps we wan't calibration data from yesterday
#TODAY=$YESTERDAY

date -u

colibri_full_redux
wait

# do some cleanup
cd $REDUX_BASE_DIR
for dir in `ls -d ${TODAY}/20*-1*/* ${TODAY}/20*-000[1,2,3]/*`; do
    cd $dir
    echo "Cleaning up $dir"
    rm -rf bias_colibri* dark_colibri* flat_colibri* 2>/dev/null
    for ps1file in `ls */catalogs/ps1_*/*fz 2>/dev/null`; do
        pfile=`basename $ps1file`
        cell=`echo $pfile | awk -F. '{print $4}'`
        scell=`echo $pfile | awk -F. '{print $5}'`
        ldir=$REDUX_BASE_DIR/ps1_stacks/$cell/$scell
        [ -d $ldir ] || mkdir -p $ldir
        [ -f "$ldir/$pfile" ] || cp $ps1file $ldir
    done
    rm -rf */biases/* */darks/* */flats/* */catalogs/* 2>/dev/null
    cd $REDUX_BASE_DIR
done

cd $REDUX_WORK_BASE_DIR
direc_str=`echo $REDUX_BASE_DIR | tr '/' '_'`
rm -r ${direc_str}_*
