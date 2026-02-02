#!/bin/sh

afile=$1
shift

fwhm=3.0
kernel_order=2
bg_order=0
sz=1024
do_phot=no
catdir=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

gflag=`echo $fwhm | awk '{ds=$1/2.355; if (ds<0.5) ds=0.5; printf("3 6 %f 4 %f 2 %f\n",0.5*ds,ds,2*ds)}'`

abase=${afile%'.fits'}
afile_wt=${abase}.wt.fits
afile_rms=${abase}.rms.fits

sat=`gethead SATURATE $afile`

[ -f $afile_rms ] || weight2rms.py $afile_wt $afile $afile_rms

ra0=`gethead ALRA $afile | awk '{printf("%.8f\n",$1)}'` ; dec0=`gethead ALDE $afile | awk '{printf("%.8f\n",$1)}'`
[ "$ra0" ] || ra0=`gethead CRVAL1 $afile | awk '{printf("%.8f\n",$1)}'`
[ "$dec0" ] || dec0=`gethead CRVAL2 $afile | awk '{printf("%.8f\n",$1)}'`
filter=`gethead FILTER $afile`

filters="g,r,i,z,y"
available=`echo $filters | grep $filter`
[ "$available" ] || { echo "That filter is not available in PS1, exiting." ; exit 1 ; }

im_url=`get_ps1imagelist.py $ra0 $dec0 $sz $filter`
wim_url=`echo $im_url | sed -e 's/\.unconv\.fits/\.unconv\.wt\.fits/g'`

bfile=ps1_${filter}.fits
bbase=${bfile%'.fits'}
bfile_var=${bbase}.var.fits
bfile_wt=${bbase}.wt.fits
bfile_rms=${bbase}.rms.fits

here=`pwd`
[ "$catdir" ] && cd $catdir
#[ -f "${bfile}.fz" ] || wget $im_url -O ${bfile}.fz 1>/dev/null 2>&1 &
#[ -f "${bfile_var}.fz" ] || wget $wim_url -O ${bfile_var}.fz 1>/dev/null 2>&1 &
[ -f "$bfile" ] || wget $im_url -O $bfile 1>/dev/null 2>&1 &
[ -f "$bfile_var" ] || wget $wim_url -O $bfile_var 1>/dev/null 2>&1 &
wait
#[ -f "${bfile}.fz" ] || { echo "Unable to download ${bfile}.fz" ; exit 1 ; }
#[ -f "${bfile_var}.fz" ] || { echo "Unable to download ${bfile_var}.fz" ; exit 2 ; }
[ -f "$bfile" ] || { echo "Unable to download $bfile" ; exit 1 ; }
[ -f "$bfile_var" ] || { echo "Unable to download $bfile_var" ; exit 2 ; }
#funpack -O ${here}/$bfile ${bfile}.fz
#funpack -O ${here}/$bfile_var ${bfile_var}.fz
if [ "$catdir" ]; then
    cp $bfile $bfile_var $here
    cd $here
fi

ps1scale.py $bfile
var2weight.py $bfile_var $bfile_wt

ls $bfile > list$$.txt
imhead -f $afile > ref_${bbase}.head
sargs="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE Y -SUBTRACT_BACK N -WEIGHT_SUFFIX .wt.fits -WEIGHT_TYPE MAP_WEIGHT -BLANK_BADPIXELS N -GAIN_DEFAULT 1.e9"
swarp @list$$.txt $sargs -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME ref_${bbase}.fits -WEIGHTOUT_NAME ref_${bbase}.wt.fits
rm ref_${bbase}.head 2>/dev/null
weight2rms.py ref_${bbase}.wt.fits ref_${bbase}.fits ref_${bbase}.rms.fits

echo "
from astropy.io.fits import getdata
from numpy import where
rms = getdata('ref_${bbase}.rms.fits')
i,j = where(rms>0)
print (\"-gd\",j.min()+1,j.max()+1,i.min()+1,i.max()+1)" | python3 > good$$.txt
gd=
[ -f "good$$.txt" ] && gd=`cat good$$.txt`

echo hotpants -inim $afile -ini $afile_rms -tni ref_$bfile_rms -tmplim ref_$bfile -outim ${abase}.diff0.fits -oni ${abase}.diff0.rms.fits -iu $sat -il -$sat -ko $kernel_order -bgo $bg_order -ng $gflag -v 0 -n i -c t $gd -oci ${bbase}.sm.fits
hotpants -inim $afile -ini $afile_rms -tni ref_$bfile_rms -tmplim ref_$bfile -outim ${abase}.diff0.fits -oni ${abase}.diff0.rms.fits -tu 1.e99 -tl -1.e99 -iu $sat -il -$sat -ko $kernel_order -bgo $bg_order -ng $gflag -v 0 -n t -c t $gd -oci ${bbase}.sm.fits 2> hotpants$$.txt

echo "
from astropy.io.fits import getdata,writeto
x,hdr=getdata('${abase}.diff0.fits',header=True)
x0=getdata('ref_$bfile')
dx=getdata('${abase}.diff0.rms.fits')
w1=getdata('$afile_wt')
w2=getdata('ref_$bfile_wt')
s=getdata('${bbase}.sm.fits')

j = s<0
x[j] += s[j]

norm = 1.*dx
j=x0>dx
norm[j] = x0[j]
j = (norm>0)*(w1>0)*(w2>0)
x[j] /= norm[j]
x[~j]=0
dx[~j]=0

writeto('${abase}.diff.rms.fits',dx,hdr,overwrite=True)
writeto('${abase}.diff.fits',x,hdr,overwrite=True)" | python3

if [ -f psf_grid.fits -a "$do_phot" = "yes" ]; then
    # need to fix normalization with respect to ps1 reference
    cp $afile_wt ${abase}.diff.wt.fits
    cp $afile_wt ${abase}.diff0.wt.fits
    [ -d ${abase}_diff_dir ] || mkdir ${abase}_diff_dir
    [ -d ${abase}_diff0_dir ] || mkdir ${abase}_diff0_dir
    cp ${abase}_dir/catalog.fits ${abase}_diff_dir
    cp ${abase}_dir/catalog.fits ${abase}_diff0_dir
    npsf=`gethead NBINX psf_grid.fits`
    x0=`gethead CRPIX1 $afile`; y0=`gethead CRPIX2 $afile`
    echo calc_psf ${abase}.diff.fits ${abase}.diff.wt.fits ${abase}_dir/seg.fits ${abase}_diff_dir/catalog.fits ${abase}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
    calc_psf ${abase}.diff.fits ${abase}.diff.wt.fits ${abase}_dir/seg.fits ${abase}_diff_dir/catalog.fits ${abase}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0 &
    echo calc_psf ${abase}.diff0.fits ${abase}.diff0.wt.fits ${abase}_dir/seg.fits ${abase}_diff0_dir/catalog.fits ${abase}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
    calc_psf ${abase}.diff0.fits ${abase}.diff0.wt.fits ${abase}_dir/seg.fits ${abase}_diff0_dir/catalog.fits ${abase}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0 &
    wait
    cat2radec.py ${abase}_diff_dir/catalog.fits ${abase}_diff_dir/${abase}_diff_radec.txt
    cat2radec.py ${abase}_diff0_dir/catalog.fits ${abase}_diff0_dir/${abase}_diff0_radec.txt

    mag0=`gethead MAG0 $afile`; vmag0=`gethead DMAG0 $afile | awk '{print $1*$1}'`
    grep -v '#' ${abase}_diff_dir/${abase}_diff_radec.txt | awk '{print $1,$2}' > rd$$.txt
    grep -v '#' ${abase}_diff_dir/${abase}_diff_radec.txt | awk '{print $7,$8,$9,$10,$11,$12,$13,$14,$15}' > rest$$.txt
    grep -v '#' ${abase}_diff_dir/${abase}_diff_radec.txt | awk '{print $3,$4,$5,$6}' > m$$.txt
    grep -v '#' ${abase}_diff0_dir/${abase}_diff0_radec.txt | awk '{print $3,$4,$5,$6}' > m0$$.txt
    paste -d ' ' m$$.txt m0$$.txt | awk '{if($1>0 && $6<0.21715) printf("%.6f %.6f %.6f %.6f\n",$5+'$mag0',sqrt($6*$6+'$vmag0'),0,0); else printf("%.6f %.6f %.6f %.6f\n",0,0,0,0)}' > m1$$.txt

    pfile=${abase}_diff_dir/${abase}_diff_radec.txt.photometry.txt
    grep '#' ${abase}_diff_dir/${abase}_diff_radec.txt > $pfile
    paste -d ' ' rd$$.txt m1$$.txt rest$$.txt >> $pfile
    rm rd$$.txt m$$.txt m0$$.txt m1$$.txt rest$$.txt
    rm ${abase}.diff.wt.fits ${abase}.diff0.wt.fits
fi

rm list$$.txt hotpants$$.txt good$$.txt
rm $afile_rms $bfile $bfile_var $bfile_wt ref_$bfile ref_$bfile_wt ref_$bfile_rms ${bbase}.sm.fits 
