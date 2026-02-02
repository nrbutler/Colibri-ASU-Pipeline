#!/bin/sh

file=$1
shift

max_off=2.0   #  maximum expected true offset from estimated image center (degrees)
pscale=2.01   # arcsec per pixel, unbinned

cat_dir=${REDUX_BASE_DIR}/catalogs/usno

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

overfac=`echo $max_off | awk '{print 0.5*$1+1.0}'`

cam=`gethead CCD_NAME $file | cut -c1-2`
[ "$cam" ] || cam=`basename $file | cut -c16-17`
year=`gethead DATE-OBS $file | cut -c1-4`
[ "$year" ] || year=2018

bin=`gethead BINNING $file`
[ "$bin" ] || bin=1

ra=`gethead CRVAL1 $file`
dec=`gethead CRVAL2 $file`

dx=`gethead NAXIS1 $file`
dy=`gethead NAXIS2 $file`

cdec=`echo $dec | awk '{print cos($1/57.296)}'`
dra=`echo $dx $pscale $overfac $cdec | awk '{print $1*$2*$3/$4/3660.}'`
ddec=`echo $dy $pscale $overfac | awk '{print $1*$2*$3/3600.}'`

ras=`echo "$ra $dra"  | awk '{r0=$1-$2/2+360; r1=$1+$2/2+360; r00 = int(r0/10)*10+5; r11 = int(r1/10)*10+5; for (i=r00;i<=r11;i+=10) if (i<360) print i; else print i-360}'`
decs=`echo "$dec $ddec"  | awk '{d0=$1-$2/2+180; d1=$1+$2/2+180; d00 = int(d0/2)*2+1; d11 = int(d1/2)*2+1; for(i=d00;i<=d11;i+=2) if(i<90) print i; else print i-180}'`

for ra0 in $ras ; do
    for dec0 in $decs ; do
        echo grab_usno_local.py ${cat_dir}/usno_b1_${ra0}_${dec0}.fits.gz usno_b1_${cam}_${ra0}_${dec0}.fits $ra $dra $dec $ddec $year
        grab_usno_local.py ${cat_dir}/usno_b1_${ra0}_${dec0}.fits.gz usno_b1_${cam}_${ra0}_${dec0}.fits $ra $dra $dec $ddec $year &
    done
done
wait

echo "
from astropy.io.fits import getdata,writeto,getheader
from glob import glob
from numpy import zeros
files=glob('usno_b1_${cam}_*_*[0-9].fits')
n=0
for file in files: n += getheader(file)['NAXIS1']
data = zeros((4,n),dtype='float32')
n0=0
for file in files:
    n = getheader(file)['NAXIS1'] + n0
    data[:,n0:n] = getdata(file)
    n0=n

h=(data[2]<16)*(data[3]<999)
writeto('usno_radec16_${cam}.fits',data[:,h])
writeto('usno_radec_${cam}.fits',data)

files=glob('usno_b1_${cam}_*_*[0-9].pm.fits')
n=0
for file in files: n += getheader(file)['NAXIS1']
if (n>0):
    data = zeros((6,n),dtype='float32')
    n0=0
    for file in files:
        n = getheader(file)['NAXIS1'] + n0
        data[:,n0:n] = getdata(file)
        n0=n

    writeto('usno_radec_${cam}.pm.fits',data)" | python3

rm usno_b1_${cam}_*.fits 2>/dev/null
