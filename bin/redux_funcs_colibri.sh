#!/bin/bash

export REDUX_BASE_DIR=/usr/local/var/colibri
export REDUX_WORK_BASE_DIR=/home/colibri
export PATH=${REDUX_BASE_DIR}/bin:/usr/local/astrometry/bin:${REDUX_BASE_DIR}/python_modules:/usr/bin:/usr/local/bin:/bin
export PYTHONPATH=${REDUX_BASE_DIR}/python_modules

export REDUX_LOCKFILE=${REDUX_BASE_DIR}/colibri.lock
export REDUX_LOCKFILE1=${REDUX_BASE_DIR}/colibri.lock1
export REDUX_LOCKFILE2=${REDUX_BASE_DIR}/colibri.lock2

export SEXTRACTOR_DIR=${REDUX_BASE_DIR}/sextractor
export SWARP_DIR=${REDUX_BASE_DIR}/swarp
export ASTNET_DIR=${REDUX_BASE_DIR}/astnet
export CALFILE_DIR=${REDUX_BASE_DIR}/calfiles

# set test=echo for testing purposes
#test=echo
test=

# all functions look at a day's worth of data defined by TODAY:
export TODAY=`date -u +20%y%m%d`
export YESTERDAY=`date -d "yesterday" -u +20%y%m%d`

# location of raw data and summary (web) files
export raw_data_archive=/nas/archive-colibri/raw
export web_pc=tcs-a
export web_user=reducer
export web_data_dir=/usr/local/var/www/main/colibri

nfile_min=2
nflat_min=6
nbias_min=6

export NBATCH=`grep processor /proc/cpuinfo | wc -l`
NBATCH=`echo $NBATCH | awk '{print $1/2}'`

cd $REDUX_BASE_DIR

function colibri_setdirs() {
    # default directories, can be changed manually
    $test rsync -a -f"+ */" -f"- *" --chmod="g=rwx" ${raw_data_archive}/${TODAY}/executor/images/ $REDUX_BASE_DIR/${TODAY} 2>/dev/null
    export BIAS_DIRS=`ls -d ${TODAY}/20*-0002/1/* 2>/dev/null`   # Unbinned only
    export DARK_DIRS=`ls -d ${TODAY}/20*-0003/1/* 2>/dev/null`   # Unbinned only
    export FLAT_DIRS=`ls -d ${TODAY}/20*-0001/*/* 2>/dev/null`
    export STANDARD_DIRS=`ls -d ${TODAY}/20*-0005/*/*/ ${TODAY}/20*-2*/*/*/ ${TODAY}/20*-0[1-9]*/*/*/  2>/dev/null`
    export GRB_DIRS=`ls -d ${TODAY}/20*-1*/*/* 2>/dev/null`
    env | grep _DIRS | grep -v XDG
}

function colibri_copy_files() {
    # populate a data directory with raw images from the data server
    for dir in $source_list; do
        today=`echo $dir | awk -F/ '{print $1}'`
        dir0=`echo $dir | sed -e "s/${today}\///g"`
        rm $dir/20*.fits.* 2>/dev/null
        $test ln -s ${raw_data_archive}/${today}/executor/images/$dir0/20*.fits.* $dir/ 2>/dev/null
    done
}

function colibri_create_manifest() {
    # group the raw fits files for reduction
    for dir in $source_list; do
        n=`ls $dir | grep fits.fz | wc -l`
        [ "$n" -eq 0 ] && continue
        $test cd $dir
        $test ls 20*[d,o].fits.fz | awk '{cam=substr($1,16,2); print $1>cam"_list.txt"}'
        $test cd $REDUX_BASE_DIR
    done
}

function colibri_create_multi_manifest() {
    # group the raw fits files for reduction
    days=`for dir in $source_list; do echo $dir | awk -F/ '{print $1}' ; done | sort -u`
    for day in $days ; do
        dirs=`for dir in $source_list; do echo $dir | grep $day ; done`
        dir0=`echo $dirs | awk '{print $1}'`
        oid=`echo $dir0 | awk -F/ '{print $2}'`
        vid=`echo $dir0 | awk -F/ '{print $3}'`
        cd ${day}/${oid}/$vid
        n=`ls */20*[f,d,b].fits.fz | wc -l`
        [ "$n" -eq 0 ] && continue
        $test ls */20*[f,d,b].fits.fz | awk '{split($1,ar,"/"); cam=substr(ar[2],16,2); print $1>cam"_list.txt"}'
        $test cd $REDUX_BASE_DIR
    done
}

function cam_redux() {
    local list=$1
    local ns=$2
    for filter in `gethead -f FILTER @$list | sort -u`; do
        echo redux_colibri $list filter=$filter cleanup=no full_stack=yes ns=$ns
        $test redux_colibri $list filter=$filter cleanup=no full_stack=yes ns=$ns &
    done
    wait
}

function colibri_do_redux() {
    local ns=$1
    [ "$ns" ] || ns=4096
    # science data reduction (after bias and flat construction)
    echo "Reducing frames: $source_list"
    colibri_copy_files
    colibri_create_manifest
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir/2*fits.fz | wc -l`
            [ "$n" -lt "$nfile_min" ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            [ -f nfiles_last.txt ] && cp nfiles_last.txt nfiles_last_redux.txt
            for list in `ls C?_list.txt 2>/dev/null`; do
                $test cam_redux $list $ns &
            done
            wait
            workdir=`cat workdir.txt`
            $test cd $workdir
            $test make_full_summary.sh savedir=$REDUX_BASE_DIR/$dir &
            $test cd $REDUX_BASE_DIR
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function colibri_do_opt_redux() {
    local cleanup=$1
    # science data reduction (after bias and flat construction)
    echo "Reducing frames: $source_list"
    colibri_copy_files
    colibri_create_manifest
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir/2*fits.fz | wc -l`
            [ "$n" -lt "$nfile_min" ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            [ -f nfiles_last.txt ] && cp nfiles_last.txt nfiles_last_redux.txt
            for list in `ls C?_list.txt 2>/dev/null`; do
                for filter in `gethead -f FILTER @$list | sort -u`; do
                    echo redux_colibri $list filter=$filter do_nir_sky=no cleanup=$cleanup
                    $test redux_colibri $list filter=$filter do_nir_sky=no cleanup=$cleanup
                done
            done
            $test cd $REDUX_BASE_DIR
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function colibri_do_bias() {
    # create bias frames and store to bias bank, work in parallel
    source_list=$BIAS_DIRS
    echo "Bias frames: $source_list"
    colibri_copy_files
    colibri_create_multi_manifest
    days=`for dir in $source_list; do echo $dir | awk -F/ '{print $1}' ; done | sort -u`
    for day in $days ; do
        dirs=`for dir in $source_list; do echo $dir | grep $day ; done`
        dir0=`echo $dirs | awk '{print $1}'`
        oid=`echo $dir0 | awk -F/ '{print $2}'`
        vid=`echo $dir0 | awk -F/ '{print $3}'`
        cd ${day}/${oid}/$vid
        for list in `ls C?_list.txt 2>/dev/null`; do
	    for rbias in `gethead -p BINNING READMODE @$list | awk '{print $2":"$3}' | sort -u`; do
		bin=`echo $rbias | awk -F: '{print $1}'`
		rmode=`echo $rbias | awk -F: '{print $2}'`
                echo bias_colibri $list bin=$bin rmode=$rmode
                $test bias_colibri $list bin=$bin rmode=$rmode &
            done
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function colibri_do_dark() {
    # create dark frames and store to dark bank, work in parallel
    source_list=$DARK_DIRS
    echo "Dark frames: $source_list"
    colibri_copy_files
    colibri_create_manifest
    days=`for dir in $source_list; do echo $dir | awk -F/ '{print $1}' ; done | sort -u`
    for day in $days ; do
        dirs=`for dir in $source_list; do echo $dir | grep $day ; done`
        for dir0 in $dirs; do
            cd $dir0
            for list in `ls C?_list.txt 2>/dev/null`; do
	         for rbias in `gethead -p BINNING READMODE @$list | awk '{print $2":"$3}' | sort -u`; do
                    bin=`echo $rbias | awk -F: '{print $1}'`
                    rmode=`echo $rbias | awk -F: '{print $2}'`
                    echo dark_colibri $list bin=$bin rmode=$rmode
                    $test dark_colibri $list bin=$bin rmode=$rmode &
                done
            done
            wait
            cd $REDUX_BASE_DIR
        done
    done
}

function colibri_do_flat() {
    # create flats frames and store to flat bank, work in parallel
    source_list=$FLAT_DIRS
    echo "Flat frames: $source_list"
    colibri_copy_files
    colibri_create_multi_manifest
    days=`for dir in $source_list; do echo $dir | awk -F/ '{print $1}' ; done | sort -u`
    for day in $days ; do
        dirs=`for dir in $source_list; do echo $dir | grep $day ; done`
        dir0=`echo $dirs | awk '{print $1}'`
        oid=`echo $dir0 | awk -F/ '{print $2}'`
        vid=`echo $dir0 | awk -F/ '{print $3}'`
        cd ${day}/${oid}/$vid
        for list in `ls C?_list.txt 2>/dev/null`; do
	    for bfilter in `gethead -p BINNING FILTER @$list | awk '{print $2":"$3}' | sort -u`; do
                bin=`echo $bfilter | awk -F: '{print $1}'`
                filter=`echo $bfilter | awk -F: '{print $2}'`
                echo flat_colibri $list bin=$bin filter=$filter
                $test flat_colibri $list bin=$bin filter=$filter &
            done
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function colibri_do_standards() {
    # set of commands to do all standard stars
    source_list=$STANDARD_DIRS
    [ "$source_list" ] && colibri_do_redux
}

function colibri_full_redux() {
    # do everything
    echo "Doing full redux for TODAY=$TODAY"
    colibri_setdirs
    colibri_do_bias
    colibri_do_dark
    colibri_do_flat
    #colibri_do_standards
}

#function alan() {
    #. /home/alan/colibri/profile
#}
