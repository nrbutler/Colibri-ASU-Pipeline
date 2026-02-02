#!/bin/sh

file=$1
shift

cat_dir=${REDUX_BASE_DIR}/catalogs/ps1

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

year=`gethead DATE-OBS $file | cut -c1-4`
[ "$year" ] || year=2020

cam=`gethead CCD_NAME $file | cut -c1-2`
[ "$cam" ] || cam=`basename $file | cut -c16-17`

# need min and max dec
nx=`gethead NAXIS1 $file`
ny=`gethead NAXIS2 $file`
echo "$nx $ny" | awk 'BEGIN{n=5}{dx=($1-1)/n;dy=($2-1)/n;for(i=0;i<=n;i++) for (j=0;j<=n;j++) print 1+i*dx,1+j*dy}' > xy$$.txt
xy2radec.py $file xy$$.txt > radec$$.txt
dec_min=`sort -n -k 2 radec$$.txt | awk '{d=$2; d0=int(d/2)*2+1; if(d<0) d0=d0-2; if(d0<-90) d0=d0+2; print d0;exit}'`
dec_max=`sort -nr -k 2 radec$$.txt | awk '{d=$2; d0=int(d/2)*2+1; if(d<0) d0=d0-2; if(d0>90) d0=d0-2; print d0;exit}'`
rm xy$$.txt radec$$.txt 2>/dev/null

# get central ra and dec
echo "$nx $ny" | awk '{print 0.5*(1+$1),0.5*(1+$2)}' > xy$$.txt
xy2radec.py $file xy$$.txt > radec$$.txt
ra0=`awk '{print $1;exit}' radec$$.txt`
dec0=`awk '{print $2;exit}' radec$$.txt`
rm xy$$.txt radec$$.txt

# get ra range at center
diag=`echo "$nx $ny" | awk '{print sqrt($1*$1+$2*$2)}'`
cdec0=`echo $dec0 | awk '{print cos($1/57.29578)}'`
dra=`gethead CD1_1 CD1_2 $file | awk '{ps=sqrt($1*$1+$2*$2); print '$diag'*ps/'$cdec0'}'`
ddec=`echo "$dec_max $dec_min" | awk '{print 2+$1-$2}'`
dec0=`echo "$dec_max $dec_min" | awk '{print 0.5*($1+$2)}'`

# the dec range to cover
decs=`echo "$dec_min $dec_max" | awk '{for(i=$1;i<=$2;i+=2) print i}'`

for dec in $decs; do
    dra0=`echo $dec | awk '{cd=cos(($1-0.9)/57.29578);cd2=cos(($1+0.9)/57.29578);if(cd2<cd) cd=cd2; print '$dra'*'$cdec0'/cd}'`
    dra1=`echo "$ra0 $dra0"  | awk '{if ($2>360) print 360; else {r0=($1-$2/2)%360; r1=($1+$2/2)%360; if (r0>r1) r0=r0-360; dra=r1-r0; if(dra>360) dra=360; print dra}}'`
    ras=`echo "$ra0 $dra1"  | awk '{r0=($1-$2/2) + 360; r1=($1+$2/2) + 360; r00 = int(r0/10)*10+5; r11 = int(r1/10)*10+5; for (i=r00;i<=r11;i+=10) print i%360}' | sort -nu`
    for ra in $ras; do
        echo grab_ps1_local.py ${cat_dir}/ps1_dr1_${ra}_${dec}.fits.gz ps1_dr1_${ra}_${dec}_${cam}.fits $ra0 $dra1 $dec0 $ddec $year
        grab_ps1_local.py ${cat_dir}/ps1_dr1_${ra}_${dec}.fits.gz ps1_dr1_${ra}_${dec}_${cam}.fits $ra0 $dra1 $dec0 $ddec $year &
    done
done
wait

echo "
from astropy.io.fits import getdata,writeto,getheader
from glob import glob
from numpy import zeros
files=glob('ps1_dr1_*_*[0-9]_${cam}.fits')
n=0
for file in files: n += getheader(file)['NAXIS1']
data = zeros((4,n),dtype='float32')
n0=0
for file in files:
    n = getheader(file)['NAXIS1'] + n0
    data[:,n0:n] = getdata(file)
    n0=n

writeto('ps1_dr1_radec_${cam}.fits',data,overwrite=True)" | python3

n=`gethead NAXIS1 ps1_dr1_radec_${cam}.fits`
if [ "$n" -eq 0 ]; then
    rm ps1_dr1_radec_${cam}.fits
else
    radec2radec.py $file ps1_dr1_radec_${cam}.fits
fi

n=`gethead NAXIS1 ps1_dr1_radec_${cam}.fits`
[ "$n" -eq 0 ] && rm ps1_dr1_radec_${cam}.fits

rm ps1_dr1_*_*_${cam}.fits 2>/dev/null
