#!/bin/bash

# read in 1 or two file names (comma separated, no spaces)
IFS=',' read -ra files <<< "$1"

file=`readlink -f ${files[0]}`
shift

if [ "${#files[@]}" -gt 1 ]; then
    file1=`readlink -f ${files[1]}`
else
    file1=$file
fi

inst=$1
if [ "$inst" ]; then
   shift
else
    inst=nircam
fi
[ "$inst" ] || inst=nircam

srcfile=`echo $1 | awk -F= '/srcfile/{print $2}'`
[ "$srcfile" ] && shift

extras=$@

base=${file%'.fits'}
base1=${file1%'.fits'}
wt_file=${base}.wt.fits
wt_file1=${base1}.wt.fits
[ -f "$wt_file" ] || wt_file=${base}.weight.fits
[ -f "$wt_file1" ] || wt_file1=${base1}.weight.fits

file0=`basename $file1`
base=${file0%'.fits'}
base_dir=`dirname $file1`
dir=${base}_dir
[ -d $dir ] || mkdir $dir

cd $dir

insta=$inst
if [ "$srcfile" ]; then
    [ -f sky.list ] || radec2xy_ddoti.py $file ${base_dir}/$srcfile > sky.list
fi
if [ -f sky.list ]; then
    echo "Running in assoc mode using sky.list file"
    insta=${inst}_assoc
fi

gain=`gethead GAIN $file`
gain1=`gethead GAIN $file1`
[ "$gain" ] || {
    gain=`gethead EGAIN $file`
    [ "$gain" ] && sethead GAIN=$gain $file
    [ "$gain" ] || gain=1.0
}
[ "$gain1" ] || {
    gain1=`gethead EGAIN $file1`
    [ "$gain1" ] && sethead GAIN=$gain1 $file1
    [ "$gain1" ] || gain=1.0
}

t0=`gethead DATE-OBS $file1 | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
t1=`gethead DATE-OBE $file1 | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
[ "$t1" ] || t1=$t0

ccd=`gethead CCD_NAME $file1`
[ "$ccd" ] || ccd=`echo $file1 | cut -c16-17`

sat_level=`gethead SATURATE $file`
[ "$sat_level" ] || sat_level=60000.0
extras="-SATUR_LEVEL $sat_level $extras "

sx=`gethead NAXIS1 $file1`
sy=`gethead NAXIS2 $file1`
x=`gethead CRPIX1 $file1`
y=`gethead CRPIX2 $file1`

sx0=`gethead NAXIS1 $file`
sy0=`gethead NAXIS2 $file`
x0=`gethead CRPIX1 $file`
y0=`gethead CRPIX2 $file`

if [ -f "$wt_file" ];then
    extras="-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_GAIN Y -WEIGHT_IMAGE $wt_file,$wt_file1 -WEIGHT_THRESH 0 -GAIN $gain $extras"
else
    extras="-WEIGHT_TYPE NONE $extras"
fi

#echo sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${insta}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw
#sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${insta}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw
echo sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file,$file1 $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${insta}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw
sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file,$file1 $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${insta}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw

exposure=`gethead EXPTIME $file1`
am=`gethead AIRMASS $file1`
[ "$am" ] && airmass=`gethead AM_COEF $file1 | awk '{if('$am'>1) print $1*('$am'-1); else print 0.0}'`
[ "$airmass" ] || airmass=0.0

fwhm_nstars=`quick_mode_ddoti.py catalog.fits`

sethead AMEXT=$airmass GAIN=$gain1 CCD=$ccd EXPTIME=$exposure T0=$t0 T1=$t1 SEXZERO=25.0 $fwhm_nstars catalog.fits -x 0

#cat2radec.py catalog.fits ${base}_radec.txt

fwhm0=`gethead FWHM $file1`
[ "$fwhm0" ] || sethead $fwhm_nstars $file1
