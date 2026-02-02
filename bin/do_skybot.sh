#!/bin/sh

file=$1

d0=`gethead DATE-OBS $file`
d1=`gethead DATE-OBE $file`
if [ "$d1" = "___" ]; then
    d1=$d0
    sethead DATE-OBE=$d1 $file
fi

t0=`date -ud $d0 +%s`
t1=`date -ud $d1 +%s`
jd=`echo $t0 $t1 | awk '{printf("%.6f\n", (0.5*($1+$2)-1498106591)/86400.+2457926.696655)}'`

dt_hrs=`echo "$t0 $t1" | awk '{printf("%.6f\n",($2-$1)/3600.)}'`

# get central ra and dec
nx=`gethead NAXIS1 $file`
ny=`gethead NAXIS2 $file`
echo "$nx $ny" | awk '{print 0.5*(1+$1),0.5*(1+$2)}' > xy$$.txt
xy2radec.py $file xy$$.txt > radec$$.txt
ra0=`awk '{print $1;exit}' radec$$.txt`
dec0=`awk '{print $2;exit}' radec$$.txt`
rm xy$$.txt radec$$.txt 2>/dev/null

nxy=`echo "$nx $ny" | awk '{print sqrt($1*$1+$2*$2)/2.}'`
ps=`gethead CD1_1 CD1_2 $file | awk '{print sqrt($1*$1+$2*$2)}'`
rad_deg=`echo $nxy $ps | awk '{printf("%.1f\n",$1*$2)}'`

echo "run_skybot($ra0,$dec0,$rad_deg,$jd,$dt_hrs)"

echo "from skybot import run_skybot
run_skybot($ra0,$dec0,$rad_deg,$jd,dt_hrs=$dt_hrs)" | python3 > mp_radec.txt
[ -s mp_radec.txt ] || rm mp_radec.txt
