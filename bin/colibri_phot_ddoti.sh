#!/bin/bash

file_list=$1
shift

catdir=

# source detection
do_lightcurves=yes
detect_thresh=1.0
snr_min=2.0
mask_sigma=6.0
npsf=12
refine=1   # refine psf positions (0 no refine, 1 refine)

#
do_dao=yes           # do psf fitting

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
[ -f stack_${cam}.wcs ] || { echo "No astrometry for stack_${cam}.fits, exiting..." ; exit 5 ; }

fwhm0=`gethead FWHM stack_${cam}.fits`
[ "$fwhm0" ] || fwhm0=2.0
match_rad=`echo $fwhm0 | awk '{if($1<2) print 2.0; else print $1}'`
phot_diam=`gethead PHOTDIAM stack_${cam}.fits`
phot_diam1=`echo $phot_diam | awk '{printf("%.1f\n",$1*3)}'`

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

if [ -f "${catdir}/gaia_radec_${cam}.fits" ]; then
    echo "Using ${catdir}/gaia_radec_${cam}.fits"
else
    filef=`readlink -f stack_${cam}.fits`
    cd $catdir
    echo grab_gaia_local.sh $filef
    grab_gaia_local.sh $filef
    cd $here
fi
ln -s ${catdir}/gaia_radec_${cam}.fits gaia_radec.fits 2>/dev/null

if [ -f "${catdir}/ps1_dr1_radec_${cam}.fits" ]; then
    echo "Using ${catdir}/ps1_dr1_radec_${cam}.fits"
else
    filef=`readlink -f stack_${cam}.fits`
    cd $catdir
    echo grab_ps1_local.sh $filef
    grab_ps1_local.sh $filef
    cd $here
fi
ln -s ${catdir}/ps1_dr1_radec_${cam}.fits ps1_dr1_radec.fits 2>/dev/null

if [ -f "${catdir}/2mass_radec_${cam}.fits" ]; then
    echo "Using ${catdir}/2mass_radec_${cam}.fits"
else
    filef=`readlink -f stack_${cam}.fits`
    cd $catdir
    echo grab_2mass_local.sh $filef
    grab_2mass_local.sh $filef
    cd $here
fi
ln -s ${catdir}/2mass_radec_${cam}.fits 2mass_radec.fits 2>/dev/null

match_file=gaia_radec.fits
[ -f "$match_file" ] || match_file=usno_radec.fits

# set calibration info
filter=`gethead FILTER stack_${cam}.fits`
[ "$filter" ] || filter=r
cfilter=${filter}_D
calfile=${catdir}/ps1_dr1_radec${cfilter}.txt
sethead CALFILT=$cfilter stack_${cam}.fits

#
# now do the real photometry
#
chi2_stack=no
if [ -f stack_${cam}_dir/use_full_stack.txt ]; then
    echo "Using source detection from chi2 stack"
    chi2_stack=yes
else
    [ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
    mkdir stack_${cam}_dir
    run_sex_ddoti.sh stack_${cam}.fits nircam srcfile=$match_file -DETECT_THRESH $detect_thresh -PHOT_APERTURES ${phot_diam},${phot_diam1} -CHECKIMAGE_TYPE SEGMENTATION,OBJECTS -CHECKIMAGE_NAME seg.fits,mask.fits -ASSOC_RADIUS $match_rad -ASSOCSELEC_TYPE ALL
    [ -f "catalog_matches.fits" ] && rm catalog_matches.fits
fi

echo sky2cal.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $match_file $calfile
sky2cal.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $match_file $calfile
apcorr0=stack_${cam}_dir/calibration.fits.apweight

cp stack_${cam}_dir/catalog.fits base_catalog.fits

x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
if [ "$do_dao" = "yes" ]; then
    nx=`gethead NAXIS1 $file0`
    [ "$nx" -lt 4096 ] && npsf=`echo $nx $npsf | awk '{printf("%.0f\n",$1*$2/4096)}'`
    echo ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0 find_psf_stars
    ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0 find_psf_stars
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
fi

# determine local calibration scaling factors (apweight)
echo ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0
ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0

# calibrate the stack
rm good_cal.fits 2>/dev/null
echo calibrate_ddoti.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/calibration.fits doplot
calibrate_ddoti.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/calibration.fits doplot > stack_${cam}_dir/stack_${cam}.report

if [ "$chi2_stack" = "no" ]; then
    #
    # now we do some source cleaning
    #
    #  we will throw out some new sources (ddoti_trimcat)
    #     clustered, particularly around bright stars, and with snr<snr_min
    #  we will flag new sources which should potentially be
    #    identified with known high proper motion stars (flag_high_props)
    #

    # check for faint sources matches
    if [ -f "ps1_dr1_radec.fits" ]; then
        echo match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits ps1_dr1_radec.fits $snr_min
        match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits ps1_dr1_radec.fits $snr_min
    fi
    if [ -f "usno_radec.fits" -a "$astrom_file" = "gaia_radec.fits" ]; then
        echo match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits usno_radec.fits $snr_min
        match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits usno_radec.fits $snr_min
    fi
    if [ -f "2mass_radec.fits" ]; then
        echo match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits 2mass_radec.fits $snr_min
        match2faint.py  stack_${cam}.fits stack_${cam}_dir/calibrated_catalog.fits 2mass_radec.fits $snr_min
    fi

    # trim out clustered detections, sources near edge, assign negative ids to new sources
    ref=stack_${cam}_dir/calibrated_catalog.fits
    echo ddoti_trimcat.py stack_${cam}.fits $ref $snr_min 0 AUTO
    ddoti_trimcat.py stack_${cam}.fits $ref $snr_min 0 AUTO

    #check for minor planets
    d1=`gethead DATE-OBE stack_${cam}.fits`
    if [ "$d1" = "___" ]; then
        d0=`gethead DATE-OBS stack_${cam}.fits`
        sethead DATE-OBE=$d0 stack_${cam}.fits
    fi
    do_skybot.sh stack_${cam}.fits
    [ -s mp_radec.txt ] || rm mp_radec.txt 2>/dev/null

    # mark stars close to high proper motion stars
    echo flag_high_props.py $ref unmatched_usno_radec.pm.fits mp_radec.txt
    flag_high_props.py $ref unmatched_usno_radec.pm.fits mp_radec.txt
fi

#
# do photometry on individual files
#
sethead PHOTDIAM=$phot_diam @$file_list

for file in `cat $file_list`; do
    base=${file%'.fits'}
    [ -d "${base}_dir" ] || mkdir ${base}_dir
    cp base_catalog.fits ${base}_dir/catalog.fits
    x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
    dx=`echo $x | awk '{printf("%.0f\n",'$x0'-$1+1)}'`
    dy=`echo $y | awk '{printf("%.0f\n",'$y0'-$1+1)}'`
    x1=`gethead NAXIS1 $file | awk '{printf("%.0f\n",'$dx'+$1-1)}'`
    y1=`gethead NAXIS2 $file | awk '{printf("%.0f\n",'$dy'+$1-1)}'`
    reg="$dx:$x1,$dy:$y1"
    if [ "$do_dao" = "yes" ]; then
        [ -f ${base}.psf_grid.fits ] && rm ${base}.psf_grid.fits
        fscale0=`gethead FSCALE0 $file | awk '{printf("%.0f\n",$1)}'`
        [ "$fscale0" -gt 10 ] && cp psf_grid.fits ${base}.psf_grid.fits  # don't estimate the psf when zero point is terrible
        echo calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine
        calc_psf $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits stack_${cam}.psf_stars.fits ${base}.psf_grid.fits $npsf $x0 $y0 $refine &
    else
        echo calc_phot $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits $x0 $y0
        calc_phot $file ${base}.wt.fits stack_${cam}_dir/seg.fits[$reg] ${base}_dir/catalog.fits $x0 $y0 &
    fi
    gonogo
done
wait; go_iter=0

# copy some info to the catalog files
gethead -a DATE-OBS DATE-OBE GAIN EXPTIME FWHM @$file_list | sed -e 's/\.fits//g' -e 's/-//g' -e 's/://g' | awk '{sub(/T/,"_",$2); sub(/T/,"_",$3); printf("sethead T0=%s T1=%s GAIN=%f EXPTIME=%f FWHM=%f %s_dir/catalog.fits -x 0\n",$2,$3,$4,$5,$6,$1)}' | sh

apcorr=stack_${cam}_dir/catalog.fits.apweight
if [ -f "$apcorr" ]; then
    for file in `cat $file_list`; do
        base=${file%'.fits'}
        echo ddoti_psffit.py $file ${base}_dir/catalog.fits $apcorr
        ddoti_psffit.py $file ${base}_dir/catalog.fits $apcorr &
        gonogo
    done
    wait; go_iter=0
fi

# make a deep mask, possibly for later
cat2radec.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/stack_${cam}_radec.txt
sources2mask.py stack_${cam}.wt.fits stack_${cam}_dir/mask.fits maskstack_${cam}.fits $mask_sigma

#
# calibrate the frames to the calibrated stack
#
ref=stack_${cam}_dir/calibrated_catalog.fits
for file in `cat $file_list`; do
    base=${file%'.fits'}
    echo calibrate_ddoti_corr.py ${base}_dir/catalog.fits $ref
    calibrate_ddoti_corr.py ${base}_dir/catalog.fits $ref > ${base}_dir/${base}.report &
    gonogo
done
wait; go_iter=0

# do photometry on the ps1 stack, and calibrate if present
# also calculate and normalize a difference image
pfilter=$filter
[ "$pfilter" = "B" ] && pfilter=g
[ "$pfilter" = "gri" ] && pfilter=r
[ "$pfilter" = "zy" ] && pfilter=z
ps1file=ps1_stack_${pfilter}.fits

if [ -f psf_grid.fits -a -f "$ps1file" ]; then
    calfile_ps1=${catdir}/ps1_dr1_radec${pfilter}.txt
    npsf=`gethead NBINX psf_grid.fits`
    x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
    ps1_base=${ps1file%'.fits'}
    [ -d "${ps1_base}_dir" ] && rm -r ${ps1_base}_dir
    mkdir ${ps1_base}_dir
    cp base_catalog.fits ${ps1_base}_dir/catalog.fits
    [ -f psf_grid0.fits ] && rm psf_grid0.fits 
    echo calc_psf n$ps1file n${ps1_base}.wt.fits stack_${cam}_dir/seg.fits ${ps1_base}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 1
    calc_psf n$ps1file n${ps1_base}.wt.fits stack_${cam}_dir/seg.fits ${ps1_base}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 1

    sethead CALFILT=$pfilter $ps1file n$ps1file
    sky2cal.py $ps1file ${ps1_base}_dir/catalog.fits $match_file $calfile_ps1
    apcorr0=${ps1_base}_dir/calibration.fits.apweight

    ddoti_psffit.py $ps1file ${ps1_base}_dir/catalog.fits $apcorr0
    calibrate_ddoti.py ${ps1_base}_dir/catalog.fits ${ps1_base}_dir/calibration.fits doplot > ${ps1_base}_dir/${ps1_base}.report
    maglim=`awk '/10-sigma limiting magnitude/{print $4;exit}' ${ps1_base}_dir/${ps1_base}.report`
    sethead MAGLIM=$maglim $ps1file
    ps1match.py stack_${cam}_dir/calibrated_catalog.fits ${ps1_base}_dir/calibrated_catalog.fits catalog_matches.fits 1.0

    psf_deconvolve psf_grid.fits psf_grid0.fits 
    echo psf_subtract stack_${cam}.fits stack_${cam}.wt.fits n$ps1file n${ps1_base}.wt.fits delta_psf_grid.fits \!stack_${cam}.diff.fits \!stack_${cam}.diff.wt.fits \!stack_${cam}.ref.fits 1
    psf_subtract stack_${cam}.fits stack_${cam}.wt.fits n$ps1file n${ps1_base}.wt.fits delta_psf_grid.fits \!stack_${cam}.diff.fits \!stack_${cam}.diff.wt.fits \!stack_${cam}.ref.fits 1
    [ -d "stack_${cam}.diff_dir" ] && rm -r stack_${cam}.diff_dir
    mkdir stack_${cam}.diff_dir
    cp base_catalog.fits stack_${cam}.diff_dir/catalog.fits
    echo calc_psf stack_${cam}.diff.fits stack_${cam}.diff.wt.fits stack_${cam}_dir/seg.fits stack_${cam}.diff_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
    calc_psf stack_${cam}.diff.fits stack_${cam}.diff.wt.fits stack_${cam}_dir/seg.fits stack_${cam}.diff_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0

    cp stack_${cam}_dir/catalog.fits.apweight stack_${cam}.diff_dir/catalog.fits.apweight
    mag0=`awk '/Magnitude Offset/{print $3;exit}' stack_${cam}_dir/stack_${cam}.report`
    echo calibrate_ddoti.py stack_${cam}.diff_dir/catalog.fits mag0=$mag0
    calibrate_ddoti.py stack_${cam}.diff_dir/catalog.fits mag0=$mag0
    rm subtraction_matches.fits 2>/dev/null
    ps1match.py stack_${cam}_dir/calibrated_catalog.fits stack_${cam}.diff_dir/calibrated_catalog.fits subtraction_matches.fits 1.0
fi
