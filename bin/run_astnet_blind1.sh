#!/bin/sh

file=$1
shift

radius=0.5  # in degrees
xyradius=2000.0 # pixels from center
inst=coatli
nstars=1000
thresh=20.0
# leave these blank for blind search:
ra=
dec=
x0=
y0=
parity=

cleanup=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f $file ] || { echo "No file $file" ; exit 1 ; }

[ "$ra" ] || ra=`gethead ETROBRA $file`
[ "$dec" ] || dec=`gethead ETROBDE $file`

backend=${ASTNET_DIR}/coatli_astnet_backend.cfg

[ "$parity" ] && parity="--parity $parity"

search=
[ "$ra" -a "$dec" ] && search="--ra $ra --dec $dec --radius $radius"

ps=`gethead CD1_1 CD1_2 $file | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600)}'`
pscale1=`echo $ps | awk '{printf("%.2f\n",$1*0.95)}'`
pscale2=`echo $ps | awk '{printf("%.2f\n",$1*1.05)}'`
scale="--scale-units arcsecperpix --scale-low $pscale1 --scale-high $pscale2"

nx=`gethead NAXIS1 $file`
ny=`gethead NAXIS2 $file`
[ "$x0" ] || x0=`gethead CRPIX1 $file`
[ "$y0" ] || y0=`gethead CRPIX2 $file`
[ "$x0" ] || x0=`echo $nx | awk '{printf("%.0f\n",$1/2.)}'`
[ "$y0" ] || y0=`echo $ny | awk '{printf("%.0f\n",$1/2.)}'`
crpix="--crpix-x $x0 --crpix-y $y0"

base=${file%'.fits'}
[ -d ${base}_dir ] || run_sex.sh $file $inst -DETECT_THRESH $thresh
#grep -v '#' ${base}_dir/${base}_radec.txt | sort -n -k 3 | awk '{if(NR<='$nstars') print $8,$9}' > xy$$.txt
#grep -v '#' ${base}_dir/${base}_radec.txt | sort -n -k 4 | awk 'BEGIN{n=0}{dx=$8-'$x0';dy=$9-'$y0';d=sqrt(dx*dx+dy*dy); if(d<'$xyradius' && $4>=5.e-4 && $4<1) {n=n+1; if(n<='$nstars') print $8,$9}}' | sort -n -k 1 > xy$$.txt
grep -v '#' ${base}_dir/${base}_radec.txt | sort -n -k 3 | awk 'BEGIN{n=0}{dx=$8-'$x0';dy=$9-'$y0';d=sqrt(dx*dx+dy*dy); if(d<'$xyradius' && $4>=5.e-4 && $4<1) {n=n+1; if(n<='$nstars') print $8,$9}}' | sort -n -k 1 > xy$$.txt

[ -f ${base}.wcs ] && rm ${base}.wcs

echo 'X E "" "" "" "" "" ""
Y E "" "" "" "" "" ""' > hf$$.txt

echo "from astropy.io.fits import tableload
x=tableload('xy$$.txt','hf$$.txt')
x.writeto('${base}.xy',overwrite=True)" | python3

echo "solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 10 -T --pnm out.pnm --new-fits none $parity"
solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 10 -T --pnm out.pnm --new-fits none $parity 2>&1

[ -f "${base}.wcs" ] && cphead CRVAL1 CRVAL2 CD1_1 CD1_2 CD2_1 CD2_2 CRPIX1 CRPIX2 ${base}.wcs $file

if [ "$cleanup" = "yes" ]; then
    rm ${base}.xy ${base}.axy xy$$.txt hf$$.txt ${base}.rdls ${base}.match ${base}.corr ${base}-indx.xyls ${base}.solved 2>/dev/null
fi
