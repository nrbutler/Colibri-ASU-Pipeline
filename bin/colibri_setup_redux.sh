#!/bin/bash

file0=$1
shift

#matching,wcs
ps=0.0001056 # platescale

cleanup=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f "$file0" ] || { echo "Cannot find first file $file0" ; exit 1 ; }

filter=`gethead FILTER $file0`
cam=`basename $file0 | cut -c16-17`
[ "$cam" = "C0" ] && cam=C1
tag=${cam}_$filter

here=`pwd`

catdir=${here}/catalogs
biasdir=${here}/biases
darkdir=${here}/darks
flatdir=${here}/flats
[ -d $catdir ] || mkdir $catdir
[ -d $biasdir ] || mkdir $biasdir
[ -d $darkdir ] || mkdir $darkdir
[ -d $flatdir ] || mkdir $flatdir

# copy the data into the current directory
bin=`gethead BINNING $file0`
[ "$bin" ] || bin=1
if [ "$bin" -gt 1 ]; then
    ps=`echo "$ps $bin" | awk '{printf("%.6f\n",$1*$2)}'`
fi

rmode=`gethead READMODE $file0`
[ "$rmode" ] || rmode=0
day=`basename $file0 | cut -c1-8`

#
# begin setup for flat fielding, bias subtraction, etc.
#  need day, cam, bin, rmode, filter
#

biasfile=bias_${cam}.fits
if [ -f ${biasdir}/$biasfile ]; then
    echo "Using biasfile ${biasdir}/$biasfile"
else
    bias=`find_bias.sh $day $cam $bin $rmode`
    cd $biasdir
    echo funpack -O $biasfile $bias
    funpack -O $biasfile $bias
    cd $here
fi
darkbase=dark_${cam}
darkfile0=`ls ${darkdir}/${darkbase}* 2>/dev/null | head -1`
if [ -f "$darkfile0" ]; then
    echo "Using darkfiles like  $darkfile0"
else
    darks=`find_dark.sh $day $cam $bin $rmode`
    cd $darkdir
    for dark in $darks; do
        dark0=`basename $dark`
        dt=`echo ${dark0%'.fits.fz'} | awk -F_ '{print $NF}'`
        darkfile0="${darkbase}_${dt}.fits"
        echo funpack -O $darkfile0 $dark
        funpack -O $darkfile0 $dark
    done
    cd $here
fi
flatfile=flat_${tag}.fits
if [ -f ${flatdir}/$flatfile ]; then
    echo "Using flatfile ${flatdir}/$flatfile"
else
    flat=`find_flat.sh $day $cam $filter $bin`
    cd $flatdir
    echo funpack -O $flatfile $flat
    funpack -O $flatfile $flat
    cd $here
fi
