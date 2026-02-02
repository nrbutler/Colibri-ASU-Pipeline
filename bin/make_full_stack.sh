#!/bin/bash

catdir=
cam=AL
detect_thresh=0.6
snr_min=2.0
snr_cut=5.0
max_sleep=120  # max waits (5s intervals) for full stack or other reduxes

npsf=12
refine=1   # refine psf positions (0 no refine, 1 refine)

thumb_size=401
thumb_edge=30
thumb_scale=1
thumb_logscale=no
thumb_invert=yes
i0=
j0=
fix_rotation=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

here0=`pwd`
me=redux_colibri_$cam
workdir=${here0}/$me
making=${here0}/making_full_stack.txt

if [ -f $making ]; then
   n=0
   while [ -f $making -a $n -lt $max_sleep ]; do
       echo sleep 5
       sleep 5
       n=`expr $n + 1`
   done
   exit 1
else
   touch $making
   go=`cat ${here0}/redux_colibri_C?_*/current_state.txt 2>/dev/null | awk '{if($1>0) {if($1>=3) s=s+1;n=n+1}}END{if (s>=n) print 1; else print 0}'`
   n0=0
   while [ $go -eq 0 -a $n0 -lt $max_sleep ]; do
       echo sleep 5 $go
       sleep 5
       go=`cat ${here0}/redux_colibri_C?_*/current_state.txt 2>/dev/null | awk '{if($1>0) {if($1>=3) s=s+1;n=n+1}}END{if (s>=n) print 1; else print 0}'`
       n0=`expr $n0 + 1`
   done
fi

[ "$catdir" ] || { echo "No catalog dir specified" ; rm $making ; exit 2 ; }

# could build a more accurate catalog:
#   f0 = [ Sum f^2*wt / Sum wt ]^(1/2)
#   df0 = 1/f0 * Sum f*df*wt / Sum wt
#   df0 = 1/f0 * 0.4*log(10) * Sum f^2*dm*wt / Sum wt
# mag0 =  -1.25*log10( ( Sum 10^(-0.8*mag) * wt )/Sum wt )
# dmag0 = ( Sum 10^(-0.8*(mag-mag0))*dmag*wt / Sum wt )
# with, wt ~ gethead MAG0 VAR0 stack_C?.fits | awk '{print 10^(0.8*$1)/$2}'

rm -r $workdir 2>/dev/null
mkdir $workdir
cd $workdir

# identify all sky-aligned stacks
ls ${here0}/redux_colibri_C?_*/stack_C?.wcs | sed -e 's/wcs/fits/g' > list.txt

[ -s list.txt ] || { echo "No stack files found, exiting" ; rm $making ; exit 3 ; }

filters="B,g,r,i,z,gri,y,zy,J,H"
filters=`gethead -a FILTER @list.txt | awk 'BEGIN{split("'$filters'",f,","); for (i=1;i<=length(f);i++) fdict[f[i]]=i}{print $2,fdict[$2]}' | sort -n -k 2 | awk '{s=s$1" "}END{print s}'`
filter=`echo $filters | sed -e 's/ //g'`

# 
# now do some stacking
#

sargs1="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLING_TYPE NEAREST -SUBTRACT_BACK N -WEIGHT_SUFFIX .wt.fits"

if [ "$fix_rotation" = "yes" ]; then
    swarp @list.txt $sargs1 -HEADER_ONLY Y -COMBINE_TYPE CHI-MEAN -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits
    gethead -a CD1_1 CD1_2 CD2_1 CD2_2 @list.txt | awk '{s1=s1+$2;s2=s2+$3;s3=s3+$4;s4=s4+$5}END{printf("%.8f %.8f %.8f %.8f\n",s1/NR,s2/NR,s3/NR,s4/NR)}' > cd$$.txt
    cd11=`awk '{print $1}' cd$$.txt`
    cd12=`awk '{print $2}' cd$$.txt`
    cd21=`awk '{print $3}' cd$$.txt`
    cd22=`awk '{print $4}' cd$$.txt`
    sethead CD1_1=$cd11 CD1_2=$cd12 CD2_1=$cd21 CD2_2=$cd22 stack_${cam}.fits
    imhead -f stack_${cam}.fits > stack_${cam}.head
    rm stack_${cam}.fits cd$$.txt 2>/dev/null
fi

echo swarp @list.txt $sargs1 -COMBINE_TYPE CHI-MEAN -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits
swarp @list.txt $sargs1 -COMBINE_TYPE CHI-MEAN -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits

[ -s stack_${cam}.fits ] || { echo "Combined stack stack_${cam}.fits not created, exiting." ; rm $making ; exit 4 ; }

gain_eff=`gethead -a GAIN VAR0 @list.txt | awk '{s=s+1/($2*sqrt($3))}END{print 0.5*NR/s}'`
#sethead CCD_NAME=$cam FILTER=$filter CALFILT=$filter VAR0=1.0 GAIN=1.E9 stack_${cam}.fits stack_${cam}.wt.fits
sethead CCD_NAME=$cam FILTER=$filter CALFILT=$filter VAR0=2.0 GAIN=$gain_eff stack_${cam}.fits stack_${cam}.wt.fits
sethead BSCALE=0.5 stack_${cam}.wt.fits

nx=`gethead NAXIS1 stack_${cam}.fits`
ny=`gethead NAXIS2 stack_${cam}.fits`

# defined in fits2jpg.py:
[ "$i0" ] || i0=`echo $nx $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
[ "$j0" ] || j0=`echo $ny $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`

# get a fwhm
[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
run_sex.sh stack_${cam}.fits coatli -DETECT_THRESH 20.0 -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME back.fits
fwhm=`quick_mode stack_${cam}_dir/test.cat n=NF min=0`
[ "$fwhm" ] || fwhm=2.0
phot_diam0=`echo $fwhm`
phot_diam=`echo $fwhm | awk '{printf("%.1f\n",1+int($1*3))}'`
match_rad=`echo $fwhm | awk '{if($1<2) print 2.0; else print $1}'`
phot_diam1=`echo $phot_diam | awk '{printf("%.1f\n",$1*3)}'`
imarith stack_${cam}.fits stack_${cam}_dir/back.fits sub tstack_${cam}.fits
mv tstack_${cam}.fits stack_${cam}.fits
sethead FWHM=$fwhm PHOTDIAM=$phot_diam stack_${cam}.fits

# find time overlaps
gethead -a DATE-OBS DATE-OBE @list.txt > times$$.txt
for stime in `awk '{print $2}' times$$.txt`; do
    stime0=`echo $stime | sed -e 's/-//g' -e 's/://g'`
    gps=`ut2gps.py $stime0`
    echo $stime $gps
done | sort -nr -k 2 > times1$$.txt
t0=`tail -1 times1$$.txt | awk '{print $1}'`
gt0=`head -1 times1$$.txt | awk '{print $2}'`
tstart=`head -1 times1$$.txt | awk '{print $1}' | sed -e 's/-//g' -e 's/://g'`
for stime in `awk '{print $3}' times$$.txt`; do
    stime0=`echo $stime | sed -e 's/-//g' -e 's/://g'`
    gps=`ut2gps.py $stime0`
    echo $stime $gps
done | sort -n -k 2 > times2$$.txt
tstop=`awk '{if($2>'$gt0') {print $1;exit}}' times2$$.txt | sed -e 's/-//g' -e 's/://g'`
t1=`tail -1 times2$$.txt | awk '{print $1}'`
rm times$$.txt times1$$.txt times2$$.txt 2>/dev/null

sethead DATE-OBS=$t0 DATE-OBE=$t1 stack_${cam}.fits stack_${cam}.wt.fits
sethead TCSTART=$tstart TCSTOP=$tstop @list.txt   # common (overlap) time interval

# try to make a similar stack for ps1, identify all sky-aligned stacks
ls ${here0}/redux_colibri_C?_*/nps1_stack_?.fits > list_ps1.txt
if [ -s "list_ps1.txt" ]; then
 
    imhead -f stack_${cam}.fits > ps1_stack_${filter}.head
    echo swarp @list_ps1.txt $sargs1 -COMBINE_TYPE CHI-MEAN -IMAGEOUT_NAME ps1_stack_${filter}.fits -WEIGHTOUT_NAME ps1_stack_${filter}.wt.fits
    swarp @list_ps1.txt $sargs1 -COMBINE_TYPE CHI-MEAN -IMAGEOUT_NAME ps1_stack_${filter}.fits -WEIGHTOUT_NAME ps1_stack_${filter}.wt.fits

    # background subtraction from sextractor
    run_sex.sh ps1_stack_${filter}.fits coatli -DETECT_THRESH 20.0 -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME back.fits
    fwhmp=`quick_mode ps1_stack_${filter}_dir/test.cat n=NF min=0`
    [ "$fwhmp" ] || fwhmp=2.0
    phot_diamp=`echo $fwhmp | awk '{printf("%.1f\n",1+int($1*3))}'`
    sethead FWHM=$fwhmp PHOTDIAM=$phot_diamp ps1_stack_${filter}.fits
    gain_effp=`gethead -a GAIN VAR0 @list_ps1.txt | awk '{s=s+1/($2*sqrt($3))}END{print 0.5*NR/s}'`
    sethead FILTER=$filter VAR0=2.0 GAIN=$gain_effp ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits
    sethead BSCALE=0.5 ps1_stack_${filter}.wt.fits
    imarith ps1_stack_${filter}.fits ps1_stack_${filter}_dir/back.fits sub tps1_stack.fits
    mv tps1_stack.fits ps1_stack_${filter}.fits
fi

#
# pull over some catalogs
#
here=`pwd`

grb_RA=`gethead CRVAL1 stack_${cam}.fits` ; grb_DEC=`gethead CRVAL2 stack_${cam}.fits`
ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600)}'`
rad=`gethead NAXIS1 NAXIS2 stack_${cam}.fits | awk '{r=sqrt($1*$1+$2*$2)*'$ps'/120; printf("%.1f\n",r)}'`

if [ -f "${catdir}/ps1_dr1_radec0.txt" ]; then
    echo "Using ${catdir}/ps1_dr1_radec0.txt"
else
    cd $catdir
    grab_ps1_dr1.sh ra=$grb_RA dec=$grb_DEC radius=$rad 
    cd $here
fi
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

# if there's an error region, map it to pixel coordinates
pos_err0=`gethead ALUN stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
if [ "$pos_err0" ]; then
    ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600.)}'`
    ra0=`gethead ALRA stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead ALDE stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
    pos_err=`echo $pos_err0 $ps | awk '{printf("%.1f\n",3600.*$1/$2)}'`
    pos_err_as=`echo $pos_err0 | awk '{printf("%.1f\n",3600.*$1)}'`
    echo $ra0 $dec0 > radec$$.txt
    radec2xy.py stack_${cam}.fits radec$$.txt noverify > xy$$.txt
    x00=`awk '{print $1;exit}' xy$$.txt`
    y00=`awk '{print $2;exit}' xy$$.txt`
    sethead SEXX0=$x00 SEXY0=$y00 SEXERR=$pos_err stack_${cam}.fits
    rm xy$$.txt radec$$.txt
fi

# 
# now do the photometry
#

[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
mkdir stack_${cam}_dir

run_sex_ddoti.sh stack_${cam}.fits nircam srcfile=$match_file -THRESH_TYPE ABSOLUTE -DETECT_THRESH $detect_thresh -PHOT_APERTURES ${phot_diam0},${phot_diam1} -CHECKIMAGE_TYPE SEGMENTATION,OBJECTS -CHECKIMAGE_NAME seg.fits,mask.fits -ASSOC_RADIUS $match_rad -ASSOCSELEC_TYPE ALL
echo catclean.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/seg.fits stack_${cam}.fits
catclean.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/seg.fits stack_${cam}.fits

[ -f "catalog_matches.fits" ] && rm catalog_matches.fits
calfile=${catdir}/ps1_dr1_radecr.txt
echo sky2cal.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $match_file $calfile
sky2cal.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $match_file $calfile
apcorr0=stack_${cam}_dir/calibration.fits.apweight

# determine local calibration scaling factors (apweight)
echo ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0 find_psf_stars
ddoti_psffit.py stack_${cam}.fits stack_${cam}_dir/catalog.fits $apcorr0 find_psf_stars

cp stack_${cam}_dir/catalog.fits base_catalog.fits
x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
[ "$nx" -lt 4096 ] && npsf=`echo $nx $npsf | awk '{printf("%.0f\n",$1*$2/4096)}'`
echo calc_psf stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 $refine
calc_psf stack_${cam}.fits stack_${cam}.wt.fits stack_${cam}_dir/seg.fits stack_${cam}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 $refine
echo chi2flux_fix.py stack_${cam}_dir/catalog.fits
chi2flux_fix.py stack_${cam}_dir/catalog.fits

# calibrate the stack
echo calibrate_ddoti.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/calibration.fits doplot
calibrate_ddoti.py stack_${cam}_dir/catalog.fits stack_${cam}_dir/calibration.fits doplot > stack_${cam}_dir/stack_${cam}.report

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

# evaluate photometry on the ps1_stack and check for matches
if [ -f ps1_stack_${filter}.fits ]; then
    cp base_catalog.fits ps1_stack_${filter}_dir/catalog.fits

    echo calc_psf ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits stack_${cam}_dir/seg.fits ps1_stack_${filter}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 $refine
    calc_psf ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits stack_${cam}_dir/seg.fits ps1_stack_${filter}_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 $refine
    echo chi2flux_fix.py ps1_stack_${filter}_dir/catalog.fits
    chi2flux_fix.py ps1_stack_${filter}_dir/catalog.fits

    sethead CALFILT=$filter ps1_stack_${filter}.fits
    sky2cal.py ps1_stack_${filter}.fits ps1_stack_${filter}_dir/catalog.fits $match_file $calfile
    apcorr0=${ps1_base}_dir/calibration.fits.apweight

    ddoti_psffit.py ps1_stack_${filter}.fits ps1_stack_${filter}_dir/catalog.fits $apcorr0
    calibrate_ddoti.py ps1_stack_${filter}_dir/catalog.fits ps1_stack_${filter}_dir/calibration.fits doplot > ps1_stack_${filter}_dir/ps1_stack_${filter}.report
    maglim=`awk '/10-sigma limiting magnitude/{print $4;exit}' ps1_stack_${filter}_dir/ps1_stack_${filter}.report`
    sethead MAGLIM=$maglim ps1_stack_${filter}.fits

    psf_deconvolve psf_grid.fits psf_grid0.fits 
    echo psf_subtract stack_${cam}.fits stack_${cam}.wt.fits ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits delta_psf_grid.fits \!stack_${cam}.diff.fits \!stack_${cam}.diff.wt.fits \!stack_${cam}.ref.fits 1
    psf_subtract stack_${cam}.fits stack_${cam}.wt.fits ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits delta_psf_grid.fits \!stack_${cam}.diff.fits \!stack_${cam}.diff.wt.fits \!stack_${cam}.ref.fits 1
    mkdir stack_${cam}.diff_dir
    cp base_catalog.fits stack_${cam}.diff_dir/catalog.fits
    echo calc_psf stack_${cam}.diff.fits stack_${cam}.diff.wt.fits stack_${cam}_dir/seg.fits stack_${cam}.diff_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
    calc_psf stack_${cam}.diff.fits stack_${cam}.diff.wt.fits stack_${cam}_dir/seg.fits stack_${cam}.diff_dir/catalog.fits stack_${cam}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
    echo chi2flux_fix.py stack_${cam}.diff_dir/catalog.fits
    chi2flux_fix.py stack_${cam}.diff_dir/catalog.fits

    cp stack_${cam}_dir/catalog.fits.apweight stack_${cam}.diff_dir/catalog.fits.apweight
    mag0=`awk '/Magnitude Offset/{print $3;exit}' stack_${cam}_dir/stack_${cam}.report`
    echo calibrate_ddoti.py stack_${cam}.diff_dir/catalog.fits mag0=$mag0
    calibrate_ddoti.py stack_${cam}.diff_dir/catalog.fits mag0=$mag0

    rm subtraction_matches.fits 2>/dev/null
    ps1match.py stack_${cam}_dir/calibrated_catalog.fits stack_${cam}.diff_dir/calibrated_catalog.fits subtraction_matches.fits 1.0

    echo calc_phot ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits stack_${cam}_dir/seg.fits ps1_stack_${filter}_dir/catalog.fits
    calc_phot ps1_stack_${filter}.fits ps1_stack_${filter}.wt.fits stack_${cam}_dir/seg.fits ps1_stack_${filter}_dir/catalog.fits
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

# flag chi2 sources present in reference image
if [ -f stack_${cam}.diff_dir/calibrated_catalog.fits ]; then
    echo ps1match_chi2.py stack_${cam}.diff_dir/calibrated_catalog.fits ps1_stack_${filter}_dir/calibrated_catalog.fits catalog_matches.fits
    ps1match_chi2.py stack_${cam}.diff_dir/calibrated_catalog.fits ps1_stack_${filter}_dir/calibrated_catalog.fits catalog_matches.fits
fi

#
# map catalog
#
echo "
from astropy.io.fits import open,getheader
from fit_wcs import ad2xy
hdu = open('stack_${cam}_dir/catalog.fits')
hdu0 = open('stack_${cam}_dir/calibrated_catalog.fits')
ids = hdu0[1].data['VECTOR_ASSOC']
hdu[1].data['VECTOR_ASSOC'] = ids
ra, dec = hdu[1].data['ALPHA_J2000'], hdu[1].data['DELTA_J2000']" > python$$.batch

for file in `cat list.txt`; do
    dir=`dirname $file`
    file0=`basename $file`
    base=${file0%'.fits'}
    cam1=`echo $base | awk -F_ '{print $2}'`
    rm -r ${dir}/${base}_dir 2>/dev/null
    mkdir ${dir}/${base}_dir
    touch ${dir}/${base}_dir/use_full_stack.txt
    #cp stack_${cam}.psf_stars.fits ${dir}/stack_${cam1}.psf_stars.fits
    cp catalog_matches.fits $dir
    [ -s mp_radec.txt ] && cp mp_radec.txt $dir
    [ -f usno_radec16.fits ] && cp -P usno_radec16.fits $dir
    [ -f usno_radec.fits ] && cp -P usno_radec.fits $dir
    [ -f gaia_radec.fits ] && cp -P gaia_radec.fits $dir
    [ -f ps1_dr1_radec.fits ] && cp -P ps1_dr1_radec.fits $dir
    [ -f 2mass_radec.fits ] && cp -P 2mass_radec.fits $dir
    echo "hdr = getheader('$file'); x, y = ad2xy(ra,dec,hdr)" >> python$$.batch
    echo "nx,ny = hdr['NAXIS1'],hdr['NAXIS2']" >> python$$.batch
    echo "hdu[1].data['X_IMAGE'], hdu[1].data['Y_IMAGE'] = x.clip(1,nx), y.clip(1,ny)" >> python$$.batch
    echo "hdu.writeto('${dir}/${base}_dir/catalog.fits',overwrite=True)" >> python$$.batch
done

cat python$$.batch | python3
rm python$$.batch

#
# map segmentation map
#
sargs0="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLING_TYPE NEAREST -SUBTRACT_BACK N -WEIGHT_TYPE NONE -FSCALASTRO_TYPE NONE"
for file in `cat list.txt`; do
    dir=`dirname $file`
    file0=`basename $file`
    base=${file0%'.fits'}
    imhead -f $file > ${dir}/${base}_dir/seg.head
    segfile=`readlink -f stack_${cam}_dir/seg.fits`
    cd ${dir}/${base}_dir
    echo swarp $segfile $sargs0 -IMAGEOUT_NAME seg.fits -WEIGHTOUT_NAME seg.wt.fits
    swarp $segfile $sargs0 -IMAGEOUT_NAME seg.fits -WEIGHTOUT_NAME seg.wt.fits &
    cd $here
    imhead -f $file > ${dir}/${base}_dir/mask.head
    maskfile=`readlink -f stack_${cam}_dir/mask.fits`
    cd ${dir}/${base}_dir
    echo swarp $maskfile $sargs0 -IMAGEOUT_NAME mask.fits -WEIGHTOUT_NAME mask.wt.fits
    swarp $maskfile $sargs0 -IMAGEOUT_NAME mask.fits -WEIGHTOUT_NAME mask.wt.fits &
    cd ..
    gethead -a DATE-OBS DATE-OBE GAIN EXPTIME FWHM $file0 | sed -e 's/\.fits//g' -e 's/-//g' -e 's/://g' | awk '{sub(/T/,"_",$2); sub(/T/,"_",$3); printf("sethead T0=%s T1=%s GAIN=%f EXPTIME=%f FWHM=%f %s_dir/catalog.fits -x 0\n",$2,$3,$4,$5,$6,$1)}' | sh
    cd $here
done
wait

for file in `cat list.txt`; do
    dir=`dirname $file`
    file0=`basename $file`
    base=${file0%'.fits'}
    echo calc_phot $file ${dir}/${base}.wt.fits ${dir}/${base}_dir/seg.fits ${dir}/${base}_dir/catalog.fits
    calc_phot $file ${dir}/${base}.wt.fits ${dir}/${base}_dir/seg.fits ${dir}/${base}_dir/catalog.fits &
done
wait

#
# make a summary
#
sethead BSCALE=1.0 stack_${cam}.wt.fits
ls stack_${cam}.fits > list0.txt
nfiles=`cat ../redux_colibri_C?_*/fmatched*_list.txt | wc -l`
ddoti_summary.sh list0.txt cam=${cam} savedir=`dirname $catdir` nfiles=$nfiles snr_cut=$snr_cut

rm $making
