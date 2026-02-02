#!/bin/bash

[ -f /home/colibri/MANUAL_MODE ] && { echo "file /home/colibri/MANUAL_MODE present, aborting" ; exit 1 ; }

[ -f /home/colibri/monitoring.txt ] && { echo "GRB monitor already running, exiting" ; exit 2 ; }

source /usr/local/var/colibri/bin/redux_funcs_colibri.sh
test=

#NBATCH=`echo $NBATCH | awk '{print $1/2}'`

[ -f "$REDUX_LOCKFILE" ] && { echo "lock file $REDUX_LOCKFILE present, exiting" ; exit 2 ; }

date -u

cd $REDUX_BASE_DIR
colibri_setdirs

nfiles_min=4

nfiles=0
function getnfiles() {
    local nfiles0=`ls 20*C?o.fits.fz 2>/dev/null | wc -l`
    local nfiles1=`ls 20*C1o.fits.fz 2>/dev/null | wc -l`
    local nfiles2=`ls 20*C2o.fits.fz 2>/dev/null | wc -l`
    local ar=($nfiles1 $nfiles2)
    local mn=$nfiles0
    for n in "${ar[@]}" ; do
        ((n < mn && n>0)) && mn=$n
    done
    nfiles=$mn
}

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

source_list=$GRB_DIRS
ready_source_list=

if [ "$source_list" ] ; then
   echo "GRB project initiated for $TODAY!"

   $test colibri_copy_files

   # find directories ready for processing
   for dir in $source_list ; do
       cd $dir
       getnfiles
       echo "Found >=$nfiles files each of the cameras."
       if [ "$nfiles" -lt "$nfiles_min" ]; then
           echo "Not enough data files yet"
           cd $REDUX_BASE_DIR
           continue
       fi
       [ -s "nfiles_last.txt" ] || echo 0 > nfiles_last.txt
       nfiles_last=`cat nfiles_last.txt`
       if [ "$nfiles" -le "$nfiles_last" ]; then
           echo "No new data files"
           cd $REDUX_BASE_DIR
           continue
       fi
       [ "$test" = "echo" ] || echo $nfiles > nfiles_last.txt
       if [ -s "nfiles_last_redux.txt" ]; then
           nfiles1=`cat nfiles_last_redux.txt | awk '{print '$nfiles'-$1}'`
           if [ "$nfiles1" -le "$nfiles_min" ]; then
               echo "Not enough new data files yet"
               cd $REDUX_BASE_DIR
               continue
           fi
       fi
       # now we need to check that we have nfiles_min more than the last time we reduced
       ready_source_list="$ready_source_list $dir"
       cd $REDUX_BASE_DIR
   done

   touch /home/colibri/monitoring.txt

   # now process the ready directories
   source_list=$ready_source_list
   [ "$source_list" ] && $test colibri_do_redux

   rm /home/colibri/monitoring.txt

fi
