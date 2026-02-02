#!/bin/bash

file_list=$1
shift

catdir=

# source detection
do_lightcurves=yes
match_file=usno_radec.fits
detect_thresh=1.0
mask_sigma=3.0
npsf=12
refine=1   # refine psf positions (0 no refine, 1 refine)

# base reduction
inst=coatli
ab_offset=0.0
plot_new_only=new_only   # just make uncatalogued source lightcurves (leave blank otherwise)

#
do_dao=yes           # do psf fitting
psfstar_calibrate=no # use psf stars for internal photometric calibration

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$catdir" ] || { echo "No catalog dir specified" ; exit 2 ; }

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 3 ; }
cam=`basename $file0 | cut -c17-18`

#
# starting work
#

[ -f stack_${cam}.fits ] || { echo "No stack file found, exiting..." ; exit 4 ; }

[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
[ -f stack_${cam}.rms.fits ] && rm stack_${cam}.rms.fits
run_astnet_blind.sh stack_${cam}.fits thresh=20.0
[ -f stack_${cam}.wcs ] || { echo "No astrometry for stack_${cam}.fits, exiting..." ; exit 5 ; }

fwhm0=`quick_mode stack_${cam}_dir/stack_${cam}_radec.txt n=7 min=0`
[ "$fwhm0" ] || fwhm0=2.0
match_rad=`echo $fwhm0 | awk '{if($1<2) print 2.0; else print $1}'`
phot_diam=`echo $fwhm0 | awk '{printf("%.1f\n",1+int($1*3))}'`
phot_diam1=`echo $phot_diam | awk '{printf("%.1f\n",$1*3)}'`
sethead PHOTDIAM=$phot_diam FWHM=$fwhm0 stack_${cam}.fits

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600)}'`
rad=`gethead NAXIS1 NAXIS2 stack_${cam}.fits | awk '{r=sqrt($1*$1+$2*$2)*'$ps'/120; printf("%.1f\n",r)}'`

here=`pwd`
if [ -f "${catdir}/usno_radec16_${cam}.fits" ]; then
    echo "Using ${catdir}/usno_radec16_${cam}.fits"
else
    filef=`readlink -f stack_${cam}.fits`
    cd $catdir
    echo grab_usno_local.sh $filef max_offset=0.1 pscale=$ps
    grab_usno_local.sh $filef max_offset=0.1 pscale=$ps
    cd $here
fi
ln -s ${catdir}/usno_radec_${cam}.fits usno_radec.fits 2>/dev/null
ln -s ${catdir}/usno_radec_${cam}.pm.fits usno_radec.pm.fits 2>/dev/null
ln -s ${catdir}/usno_radec16_${cam}.fits usno_radec16.fits 2>/dev/null

grb_RA=`gethead CRVAL1 stack_${cam}.fits`
grb_DEC=`gethead CRVAL2 stack_${cam}.fits`
if [ -f "${catdir}/ps1_dr1_radec.txt" ]; then
    echo "PS1 calibration files already exist."
else
    cd $catdir
    grab_ps1_dr1.sh ra=$grb_RA dec=$grb_DEC radius=$rad
    cd $here
fi

#
# possibly tune up the wcs
#
filter=`gethead FILTER stack_${cam}.fits`
[ "$filter" ] || filter=r
cfilter=${filter}_D
calfile=${catdir}/ps1_dr1_radec${cfilter}.txt
tunefile=${catdir}/ps1_dr1_radec0.txt
[ -f "$tunefile" ] || tunefile=$calfile

# possibly tune up the wcs to improve upon astrometry.net
if [ -f "$tunefile" ]; then
   run_sex.sh stack_${cam}.fits $inst -DETECT_THRESH 10.0
   echo tune_wcs.py stack_${cam}.fits stack_${cam}_dir/stack_${cam}_radec.txt $tunefile
   tune_wcs.py stack_${cam}.fits stack_${cam}_dir/stack_${cam}_radec.txt $tunefile
fi

#
# now do the real photometry
#

[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
mkdir stack_${cam}_dir
[ -d med_stack_${cam}_dir ] && rm -r med_stack_${cam}_dir
mkdir med_stack_${cam}_dir
radec2xy.py stack_${cam}.fits $match_file > stack_${cam}_dir/sky.list

run_sex_ddoti.sh stack_${cam}.fits nircam -DETECT_THRESH $detect_thresh -PHOT_APERTURES ${phot_diam},${phot_diam1} -CHECKIMAGE_TYPE SEGMENTATION,OBJECTS -CHECKIMAGE_NAME seg.fits,mask.fits -ASSOC_RADIUS $match_rad -ASSOCSELEC_TYPE ALL

#cp stack_${cam}_dir/catalog.fits base_catalog.fits

x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
if [ "$do_dao" = "yes" ]; then
    nx=`gethead NAXIS1 $file0`
    [ "$nx" -lt 4096 ] && npsf=`echo $nx $npsf | awk '{printf("%.0f\n",$1*$2/4096)}'`
    echo ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits noprior.txt find_psf_stars
    ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits noprior.txt find_psf_stars
    npsf_stars=`gethead NAXIS1 stack_${cam}.psf_stars.fits`
    sqrt_npsf_stars=`echo $npsf_stars | awk '{printf("%.0f\n",sqrt($1))}'`
    [ "$sqrt_npsf_stars" -lt "$npsf" ] && npsf=$sqrt_npsf_stars
    [ "$npsf" -lt 1 ] && npsf=1
    rm psf_grid.fits 2>/dev/null
    echo calc_psf stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 $refine
    calc_psf stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 $refine
else
    echo calc_phot stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits
    calc_phot stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits
    #cp stack_${cam}_dir/catalog.fits med_stack_${cam}_dir/catalog.fits
    #echo calc_phot med_stack_${cam}.fits med_stack_${cam}.wt.fits stack_${cam}_dir/seg.fits med_stack_${cam}_dir/catalog.fits
    #calc_phot med_stack_${cam}.fits med_stack_${cam}.wt.fits stack_${cam}_dir/seg.fits med_stack_${cam}_dir/catalog.fits
fi

# determine local calibration scaling factors (apweight)
echo ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits noprior.txt
ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits noprior.txt

#
# do photometry on individual files
#
sortcat.py stack_${cam}_dir/catalog.fits # sort catalog by flux and assign (negative) numbers to new sources
cp $file_list r$file_list
ls stack_20*[0-9].fits >> r$file_list 2>/dev/null
ls stack_${cam}a.fits stack_${cam}b.fits >> r$file_list
sethead PHOTDIAM=$phot_diam @r$file_list

for file in `cat r$file_list`; do
    base=${file%'.fits'}
    [ -d "${base}_dir" ] || mkdir ${base}_dir
    cp stack_${cam}_dir/catalog.fits ${base}_dir/catalog.fits
    x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
    dx=`echo $x | awk '{printf("%.0f\n",'$x0'-$1+1)}'`
    dy=`echo $y | awk '{printf("%.0f\n",'$y0'-$1+1)}'`
    x1=`gethead NAXIS1 $file | awk '{printf("%.0f\n",'$dx'+$1-1)}'`
    y1=`gethead NAXIS2 $file | awk '{printf("%.0f\n",'$dy'+$1-1)}'`
    reg="$dx:$x1,$dy:$y1"
    if [ "$do_dao" = "yes" ]; then
        [ -f ${base}.psf_grid.fits ] && rm ${base}.psf_grid.fits
        echo calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine
        calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine &
        #echo calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 0
        #calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 0 &
        # use the stack psf as prior:
        #echo calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine psf_grid.fits
        #calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine psf_grid.fits &
    else
        echo calc_phot $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits $x0 $y0
        calc_phot $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits $x0 $y0 &
    fi
    gonogo
done
wait; go_iter=0

apcorr=stack_${cam}_dir/catalog.fits.apweight
if [ -f "$apcorr" ]; then
    for file in `cat r$file_list`; do
        base=${file%'.fits'}
        echo ddoti_psffit.py $file ${base}_dir/catalog.fits $apcorr
        ddoti_psffit.py $file ${base}_dir/catalog.fits $apcorr &
        gonogo
    done
    wait; go_iter=0
fi

#
# now do some calibration
#

cat2radec.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/stack_${cam}_radec.txt

# make a deep mask, possibly for later
sources2mask.py stack_${cam}.wt.fits stack_${cam}_dir/mask.fits maskstack_${cam}.fits $mask_sigma

if [ -f "$calfile" ]; then
    echo do_match1.py stack_${cam}_dir/stack_${cam}_radec.txt $calfile 18.0 $ab_offset do_plot 
    do_match1.py stack_${cam}_dir/stack_${cam}_radec.txt $calfile 18.0 $ab_offset do_plot > stack_${cam}_dir/stack_${cam}.report
else
    filter=$filter0
    echo calibrate.py stack_${cam}_dir/stack_${cam}_radec.txt stack_${cam}_dir/sky.list 18.0 $ab_offset 
    calibrate.py stack_${cam}_dir/stack_${cam}_radec.txt stack_${cam}_dir/sky.list 18.0 $ab_offset > stack_${cam}_dir/stack_${cam}.report
fi

mag0=`grep 'Magnitude Offset' stack_${cam}_dir/stack_${cam}.report | awk '{print $3;exit}'`
dmag0=`grep 'Magnitude Offset' stack_${cam}_dir/stack_${cam}.report | awk '{print $5;exit}'`
sethead CALERR=$dmag0 stack_${cam}.fits stack_${cam}.wt.fits

#
# do calibration for individual frames
#
ref=stack_${cam}_radec0.xy
nx=`gethead NAXIS1 $file0`; ny=`gethead NAXIS2 $file0`
nx0=`gethead NAXIS1 stack_${cam}.fits`; ny0=`gethead NAXIS2 stack_${cam}.fits`
ifile0=stack_${cam}_dir/stack_${cam}_radecpsf.txt
if [ -f "$ifile0" -a "$psfstar_calibrate" = "yes" ]; then
    ifile=stack_${cam}_dir/stack_${cam}_radecpsf.txt.photometry.txt
    grep '#' $ifile0 > $ifile
    grep -v '#' $ifile0 | awk '{print $1,$2,$3+'$mag0',$4,$5+'$mag0',$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' >> $ifile
    grep -v '#' $ifile | awk '{print $8,$9,$NF,$3,$4}' > $ref
else
    ifile=stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt
    grep -v '#' $ifile | awk '{if($NF>0) {dx=('$nx0'-'$nx')/2.+10;dy=('$ny0'-'$ny')/2.+10;x=$8;y=$9;if(x>dx && x<'$nx0'-dx && y>dy && y<'$ny0'-dy) print x,y,$NF,$3,$4}}' > $ref
fi

function docal() {
    local file=$1
    local base=${file%'.fits'}
    local exposure=`gethead EXPTIME $file`
    local gain1=`gethead GAIN $file`
    local fwhm=`gethead FWHM $file`
    local t0=`gethead DATE-OBS $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
    local t1=`gethead DATE-OBE $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
    sethead GAIN=$gain1 EXPTIME=$exposure T0=$t0 T1=$t1 FWHM=$fwhm ${base}_dir/catalog.fits -x 0
    cat2radec.py ${base}_dir/catalog.fits ${base}_dir/${base}_radec.txt
    echo calibrate.py ${base}_dir/${base}_radec.txt $ref 18.0 0.0
    calibrate.py ${base}_dir/${base}_radec.txt $ref 18.0 0.0 > ${base}_dir/${base}.report
}

for file in `cat r$file_list`; do
    docal $file &
    gonogo
done
wait; go_iter=0

#
# record some statistics
#
lms=`egrep 'FWHM|10-sigma limiting|Zero' stack_${cam}_dir/stack_${cam}.report | awk '{printf("%f ",$4)}END{printf("\n")}' | awk '{printf("FWHM=%.2f MAGZERO=%.2f MAGLIM=%.2f\n",$1,$2,$3)}'`
sethead $lms stack_${cam}.fits stack_${cam}.wt.fits
mv stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt_dm.jpg calplot.jpg

# remember which filter we actually used to calibrate
sethead CALFILT=$cfilter stack_${cam}.fits

# divide source list into old and new (sfile and nfile)
sfile=stack_${cam}_radec0.txt
sfile_raw=stack_${cam}_radec0_raw.txt
grep '#' stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt > $sfile
grep '#' stack_${cam}_dir/stack_${cam}_radec.txt > $sfile_raw
grep -v '#' stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt | awk '{if ($3>0 && $NF>0) print}' | sort -n -k 3 >> $sfile
grep -v '#' stack_${cam}_dir/stack_${cam}_radec.txt | awk '{if ($3>0 && $NF>0) print}' | sort -n -k 3 >> $sfile_raw
nfile=stack_${cam}_radec.txt
nfile_raw=stack_${cam}_radec_raw.txt
grep '#' stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt > $nfile
grep '#' stack_${cam}_dir/stack_${cam}_radec.txt > $nfile_raw
grep -v '#' stack_${cam}_dir/stack_${cam}_radec.txt.photometry.txt | awk '{if ($3>0 && $NF<0) print}' | sort -n -k 3 >> $nfile
grep -v '#' stack_${cam}_dir/stack_${cam}_radec.txt | awk '{if ($3>0 && $NF<0) print}' | sort -n -k 3 >> $nfile_raw

#
# do some more source checking
#

# check for GRB information
file1=`sort -nr $file_list | head -1`
pos_err0=`gethead ALUN $file1 | awk '{printf("%.6f\n",$1)}'`
if [ "$pos_err0" ]; then
    ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600.)}'`
    ra0=`gethead -c ALRA $file1 | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead -c ALDE $file1 | awk '{printf("%.6f\n",$1)}'`
    pos_err=`echo $pos_err0 $ps | awk '{printf("%.1f\n",3600.*$1/$2)}'`
    echo $ra0 $dec0 > radec$$.txt
    radec2xy.py stack_${cam}.fits radec$$.txt > xy$$.txt
    x00=`awk '{print $1;exit}' xy$$.txt`
    y00=`awk '{print $2;exit}' xy$$.txt`
    sethead SEXX0=$x00 SEXY0=$y00 SEXERR=$pos_err stack_${cam}.fits
    rm xy$$.txt radec$$.txt
fi

# check against ps1
ps1file=${catdir}/ps1_dr1_radec0.txt
if [ -f "$ps1file" ]; then
    rm ${sfile}.matched.txt ${nfile}.matched.txt 2>/dev/null
    do_match1.py $sfile $ps1file > /dev/null 2>&1
    do_match1.py $nfile $ps1file > /dev/null 2>&1
fi

# record calibration data
echo stack_${cam}.fits > zlf_list$$.txt
cat r$file_list >> zlf_list$$.txt
for file in `cat zlf_list$$.txt`; do
    base=${file%'.fits'}
    rfile=${base}_dir/${base}.report
    dat=`awk '{if ($0~/Median Sextractor FWHM/) f=$4; if ($0~/Median Zero Point/) zp=$4; if($0~/Magnitude Offset/) {m0=$3;dm0=$5}; if ($0~/10-sigma limiting magnitude/) {print t1,t2,f,zp,$4,m0,dm0}}' $rfile`
    echo $dat | awk '{printf("sethead MAGZERO=%.2f MAGLIM=%.2f MAG0=%.6f DMAG0=%.6f '${base}.fits'\n",$2,$3,$4,$5)}' | sh &
    gonogo
done
wait; go_iter=0

# now track the new stars in the individual files
echo "# id ra dec mag dmag fwhm expos" > master_phot_${cam}.txt
grep -v '#' $nfile | awk '{printf("%6d %12.8f %12.8f %12.6f %12.6f %8.2f %12.6f\n",$NF,$1,$2,$3,$4,$7,$14)}' > master_phot_${cam}.tmp
grep -v '#' $sfile | awk '{printf("%6d %12.8f %12.8f %12.6f %12.6f %8.2f %12.6f\n",$NF,$1,$2,$3,$4,$7,$14)}' >> master_phot_${cam}.tmp
sort -n -k 4 master_phot_${cam}.tmp >> master_phot_${cam}.txt

# summary of images used
gethead -a IMID IMNUM IMPID DATE-OBS DATE-OBE EXPTIME FWHM MAGZERO MAGLIM MAG0 DMAG0 stack_${cam}.fits | sed -e 's/-//g' -e 's/://g' > times_$file_list
gethead -a IMID IMNUM IMPID DATE-OBS DATE-OBE EXPTIME FWHM MAGZERO MAGLIM MAG0 DMAG0 @r$file_list | sort -n -k 2 | sed -e 's/-//g' -e 's/://g' >> times_$file_list

# use these to make a summary fits file photometry.fits
trig_time=`gethead ALEVT $file0 | sed -e 's/-//g' -e 's/://g'`
[ "$trig_time" ] || trig_time=`gethead ALALT $file0 | sed -e 's/-//g' -e 's/://g'`
[ "$trig_time" ] || trig_time=`gethead DATE-OBS stack_${cam}.fits | sed -e 's/-//g' -e 's/://g'`
echo phot2fits.py times_$file_list master_phot_${cam}.txt $trig_time $cfilter
phot2fits.py times_$file_list master_phot_${cam}.txt $trig_time $cfilter

# make a data quality plot
zlf_plot.py photometry.fits
mv photometry.fits.jpg zlf_plot.jpg

if [ "$do_lightcurves" = "yes" ]; then

    #[ "$pos_err0" ] && plot_new_only=
    [ -d "lc_dir" ] && rm -r lc_dir
    mkdir lc_dir
    echo coatli_lc_plots.py photometry.fits $plot_new_only
    coatli_lc_plots.py photometry.fits $plot_new_only
fi

rm -r r$file_list zlf_list$$.txt master_phot_${cam}.tmp 2>/dev/null
