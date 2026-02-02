#!/bin/bash

file_list=$1
shift

inst=coatli
gain=6.2
detect_thresh=1.0
ab_offset=0.0
imid=0

do_phot=yes

wcs_corr=yes

#outfile name override
ofile=
maskfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" 1>&2 ; exit 1 ; }

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" 1>&2 ; exit 2 ; }
base0=${file0%'.fits'}
file0_wt=${base0}.wt.fits
[ -f $file0_wt ] || { echo "Cannot find first weight file $file0_wt" 1>&2 ; exit 3 ; }

slist=s_$$_$file_list
sort -n $file_list > $slist
file1=`head -1 $slist` ; sdte=`gethead DATE-OBS $file1`
file2=`tail -1 $slist` ; edte=`gethead DATE-OBE $file2`

file1b=`echo $sdte | sed -e 's/-//g' -e 's/://g' | awk -F. '{print $1}'`
file2b=`echo $edte | sed -e 's/-//g' -e 's/://g' | awk -F. '{print $1}'`

if [ "$ofile" ]; then
    base=${ofile%'.fits'}
else
    base=stack_${file1b}_${file2b}
    ofile=${base}.fits
fi
ofile_wt=${base}.wt.fits

cam=`gethead CCD_NAME $file0`
[ "$maskfile" ] || maskfile=maskstack0_${cam}.fits

if [ -f ${base0}.resamp.fits ]; then
    rs=
else
    rs="-RESAMPLE_SUFFIX .fits"
fi

sargs3="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE N -WEIGHT_SUFFIX .wt.fits $rs"

[ "$wcs_corr" = "yes" ] && cp stack_${cam}.head ${base}.head 2>/dev/null

n=`cat $slist | wc -l`

#echo swarp @$slist $sargs3 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME $ofile -WEIGHTOUT_NAME $ofile_wt 1>&2 2>/dev/null
#swarp @$slist $sargs3 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME $ofile -WEIGHTOUT_NAME $ofile_wt 2>/dev/null
echo swarp @$slist $sargs3 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME $ofile -WEIGHTOUT_NAME $ofile_wt > ${ofile}.slog
swarp @$slist $sargs3 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME $ofile -WEIGHTOUT_NAME $ofile_wt 2>> ${ofile}.slog

x0=`gethead CRPIX1 $ofile`; y0=`gethead CRPIX2 $ofile`
for file in `cat $slist`; do
    x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
    xmin=`gethead XMIN $file`; xmax=`gethead XMAX $file`
    ymin=`gethead YMIN $file`; ymax=`gethead YMAX $file`
    echo "$xmin $xmax $x $ymin $ymax $y"
done > xydata_$$.txt
xmin=`awk 'BEGIN{mn=1.e56}{if($1-$3<mn) mn=$1-$3}END{print mn+'$x0'}' xydata_$$.txt`
xmax=`awk 'BEGIN{mx=-1.e56}{if($2-$3>mx) mx=$2-$3}END{print mx+'$x0'}' xydata_$$.txt`
ymin=`awk 'BEGIN{mn=1.e56}{if($4-$6<mn) mn=$4-$6}END{print mn+'$y0'}' xydata_$$.txt`
ymax=`awk 'BEGIN{mx=-1.e56}{if($5-$6>mx) mx=$5-$6}END{print mx+'$y0'}' xydata_$$.txt`
rm xydata_$$.txt

x0=`gethead CRPIX1 stack0_${cam}.fits`; y0=`gethead CRPIX2 stack0_${cam}.fits`
nx0=`gethead NAXIS1 stack0_${cam}.fits`; ny0=`gethead NAXIS2 stack0_${cam}.fits`
nx=`gethead NAXIS1 $ofile`; ny=`gethead NAXIS2 $ofile`
xr=`gethead CRPIX1 $ofile | awk '{x=int($1-('$x0'));if ('$nx'-x>'$nx0') x='$nx'-'$nx0';printf("%.0f-%.0f\n",1-x,'$nx'-x)}'`
yr=`gethead CRPIX2 $ofile | awk '{y=int($1-('$y0'));if ('$ny'-y>'$ny0') y='$ny'-'$ny0';printf("%.0f-%.0f\n",1-y,'$ny'-y)}'`
getfits $maskfile $xr $yr -o ${base}.wmap_mask.fits 1>&2

sat_level=`gethead -a SATURATE @$slist | awk 'BEGIN{s=1.e56}{if($2<s) s=$2}END{print s}'`
exptime=`gethead -a EXPTIME @$slist | awk '{s=s+$2}END{printf("%f\n", s)}'`
fwhm=`gethead -a FWHM @$slist | awk '{s=s+$2}END{printf("%.2f\n", s/NR)}'`
#igain=`echo $gain $exptime | awk '{print $1*$2}'`
wgain=`gethead -a GAIN VAR0 @$slist | awk '{s=s+1/($2*$3*$3);s0=s0+1/$3}END{print s0*s0/s}'`
var0=`gethead -a VAR0 @$slist | awk '{s=s+1/$2}END{printf("%.6f\n",1/s)}'`
sky0=`gethead -a SKYLEV @$slist | awk '{n=n+1;w0=w0+1/$2}END{print n/w0}'`
sethead XMIN=$xmin YMIN=$ymin XMAX=$xmax YMAX=$ymax VAR0=$var0 SKYLEV=$sky0 SATURATE=$sat_level DATE-OBS=$sdte FWHM=$fwhm DATE-OBE=$edte EXPTIME=$exptime GAIN=$wgain $ofile

rm $slist

if [ "$do_phot" = "yes" ]; then

    [ "$wcs_corr" = "yes" ] && cphead CRVAL1 CRVAL2 CRPIX1 CRPIX2 CD1_1 CD1_2 CD2_1 CD2_2 stack_${cam}.fits $ofile

    [ -d "${base}_dir" ] && rm -r ${base}_dir
    mkdir ${base}_dir
    x=`gethead CRPIX1 $ofile`; y=`gethead CRPIX2 $ofile`
    awk '{print $1+'$x',$2+'$y',$3,$4,$5}' match_file.xy > ${base}_dir/sky.list

    run_sex.sh ${base}.fits $inst -DETECT_THRESH $detect_thresh > /dev/null 2>&1

    ref=stack_${cam}_radec0.xy
    calibrate.py ${base}_dir/${base}_radec.txt $ref 18.0 $ab_offset > ${base}_dir/${base}.report

    lms=`egrep 'FWHM|10-sigma limiting|Zero' ${base}_dir/${base}.report | awk '{printf("%f ",$4)}END{printf("\n")}' | awk '{printf("FWHM=%.2f MAGZERO=%.2f MAGLIM=%.2f\n",$1,$2,$3)}'`
    sethead $lms $ofile $ofile_wt

    mv ${base}_dir/${base}_radec.txt ${base}_radec_raw.txt
    mv ${base}_dir/${base}_radec.txt.photometry.txt ${base}_radec.txt

    rm -r ${base}_dir

fi

sethead IMNUM=$n IMID=$imid $ofile $ofile_wt

echo $ofile

rm ${base}.head 2>/dev/null
