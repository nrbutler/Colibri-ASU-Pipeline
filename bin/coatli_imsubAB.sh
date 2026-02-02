#!/bin/sh

afile=$1
shift
bfile=$1
shift

fwhm=3.0
kernel_order=2
bg_order=0

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

gflag=`echo $fwhm | awk '{ds=$1/2.355; if (ds<0.5) ds=0.5; printf("3 6 %f 4 %f 2 %f\n",0.5*ds,ds,2*ds)}'`

abase=${afile%'.fits'}
afile_wt=${abase}.wt.fits
afile_rms=${abase}.rms.fits
bbase=${bfile%'.fits'}
bfile_wt=${bbase}.wt.fits
bfile_rms=${bbase}.rms.fits

sat1=`gethead SATURATE $afile`
sat2=`gethead SATURATE $bfile`

[ -f $afile_rms ] || weight2rms.py $afile_wt $afile $afile_rms
[ -f $bfile_rms ] || weight2rms.py $bfile_wt $bfile $bfile_rms

echo hotpants -inim $afile -ini $afile_rms -tni $bfile_rms -tmplim $bfile -outim ${abase}.diff0.fits -oni ${abase}.diff0.rms.fits -tu $sat2 -iu $sat1 -tl -$sat2 -il -$sat1 -ko $kernel_order -bgo $bg_order -ng $gflag -v 0 -n t -c t -oci ${bbase}.sm.fits
hotpants -inim $afile -ini $afile_rms -tni $bfile_rms -tmplim $bfile -outim ${abase}.diff0.fits -oni ${abase}.diff0.rms.fits -tu $sat2 -iu $sat1 -tl -$sat2 -il -$sat1 -ko $kernel_order -bgo $bg_order -ng $gflag -v 0 -n t -c t -oci ${bbase}.sm.fits 2> hotpants$$.txt

echo "
from astropy.io.fits import getdata,writeto
x,hdr=getdata('${abase}.diff0.fits',header=True)
x0=getdata('$bfile')
dx=getdata('${abase}.diff0.rms.fits')
w1=getdata('$afile_wt')
w2=getdata('$bfile_wt')
s=getdata('${bbase}.sm.fits')

j = s<0
x[j] += s[j]

norm = 1.*dx
j=x0>dx
norm[j] = x0[j]
j = (norm>0)*(w1>0)*(w2>0)
x[j] /= norm[j]
x[~j]=0

writeto('${abase}.diff.fits',x,hdr,overwrite=True)" | python3

rm hotpants$$.txt $afile_rms $bfile_rms
