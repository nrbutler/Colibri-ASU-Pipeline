#!/bin/bash

file_list=$1
shift

savedir=`pwd`

# summary
thumb_size=401
thumb_edge=30
thumb_scale=1
#thumb_logscale=yes
thumb_logscale=no
thumb_invert=yes
i0=
j0=

make_fits_thumbs=no

snr_cut=10

cam=C0

t00=`date +%s`

nfiles=0

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

file0=`sort $file_list | head -1`
[ "$filter" ] || filter=`gethead FILTER $file0`
[ "$filter" ] || filter=r

if [ "$file0" = "stack_${cam}.fits" ]; then
    dte1=`gethead DATE-OBS $file0 | sed -e 's/-//g' -e 's/://g' | awk -F. '{print $1}'`
    dte2=`gethead DATE-OBE $file0 | sed -e 's/-//g' -e 's/://g' | awk -F. '{print $1}'`
    file1=$file0
    rm $file_list
    touch $file_list
else
    dte1=`echo $file0 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`
    file1=`sort $file_list | tail -1`
    dte2=`echo $file1 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`
fi

tag=${cam}_$filter
grb_RA=`gethead CRVAL1 stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
grb_DEC=`gethead CRVAL2 stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`

pdir0=photometry_${dte1}_${tag}_${nfiles}

#
# now start working, make lightcurves, make postage stamps
#  then copy data out to a webpage and supporting directories
#

#
# record calibration data
#
echo stack_${cam}.fits > zlf_list$$.txt
cat $file_list >> zlf_list$$.txt
for file in `cat zlf_list$$.txt`; do
    base=${file%'.fits'}
    rfile=${base}_dir/${base}.report
    dat=`awk '{if ($0~/Median Zero Point/) zp=$4; if($0~/Magnitude Offset/) {m0=$3;dm0=$5}; if ($0~/10-sigma limiting magnitude/) {print t1,t2,zp,$4,m0,dm0}}' $rfile`
    echo $dat | awk '{printf("sethead MAGZERO=%.2f MAGLIM=%.2f MAG0=%.6f DMAG0=%.6f '${base}.fits'\n",$1,$2,$3,$4,$5)}' | sh &
    gonogo
done
wait; go_iter=0

[ "$cam" = "AL" ] && calscatt=`gethead CALSCAT catalog_matches.fits`  # use the measurement from image subtraction if available
[ "$calscatt" ] || calscatt=`awk '/Calibration Scatter/{print $3}' stack_${cam}_dir/stack_${cam}.report`

if [ -s $file_list ]; then
    sort $file_list > s$file_list
    gethead -a DATE-OBS DATE-OBE @s$file_list > times_${file_list}.tmpa0
    t1=`head -1 times_${file_list}.tmpa0 | awk '{print $2}'`
    t2=`tail -1 times_${file_list}.tmpa0 | awk '{print $3}'`
    sethead DATE-OBS=$t1 DATE-OBE=$t2 stack_${cam}.fits
    sed -e 's/-//g' -e 's/://g' times_${file_list}.tmpa0 > times_${file_list}.tmpa
    rm times_${file_list}.tmpa0

    gethead -af EXPTIME CRPIX1 CRPIX2 MAGZERO MAGLIM MAG0 DMAG0 AIRMASS SKYLEV VAR0 FWHM @s$file_list | sed -e 's/___/0.0/g' > times_${file_list}.tmpb
    paste times_${file_list}.tmpa times_${file_list}.tmpb | sort -n -k 2 > times_${file_list}.tmp
fi

gethead -a DATE-OBS DATE-OBE stack_${cam}.fits | sed -e 's/-//g' -e 's/://g' > times_${file_list}.a
gethead -af EXPTIME CRPIX1 CRPIX2 MAGZERO MAGLIM MAG0 DMAG0 AIRMASS SKYLEV VAR0 FWHM stack_${cam}.fits | sed -e 's/___/0.0/g' > times_${file_list}.b

paste times_${file_list}.a times_${file_list}.b > times_$file_list
[ -s times_${file_list}.tmp ] && cat times_${file_list}.tmp >> times_$file_list
rm times_${file_list}.a times_${file_list}.b times_${file_list}.tmp times_${file_list}.tmpa times_${file_list}.tmpb s$file_list 2>/dev/null

# check for Swift information
[ -s $file_list ] && file1=`sort -nr $file_list | head -1`
pos_err0=`gethead ALUN $file1 | awk '{printf("%.6f\n",$1)}'`
ref=stack_${cam}_dir/calibrated_catalog.fits
if [ "$pos_err0" ]; then
    ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600.)}'`
    ra0=`gethead ALRA $file1 | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead ALDE $file1 | awk '{printf("%.6f\n",$1)}'`
    pos_err=`echo $pos_err0 $ps | awk '{printf("%.1f\n",3600.*$1/$2)}'`
    pos_err_as=`echo $pos_err0 | awk '{printf("%.1f\n",3600.*$1)}'`
    echo $ra0 $dec0 > radec$$.txt
    radec2xy.py stack_${cam}.fits radec$$.txt noverify > xy$$.txt
    x00=`awk '{print $1;exit}' xy$$.txt`
    y00=`awk '{print $2;exit}' xy$$.txt`
    sethead SEXX0=$x00 SEXY0=$y00 SEXERR=$pos_err SEXALUN=$pos_err0 SEXALRA=$ra0 SEXALDE=$dec0 stack_${cam}.fits
    sethead SEXX0=$x00 SEXY0=$y00 PS=$ps $ref -x 0
    rm xy$$.txt radec$$.txt
fi

# use the info to make a summary fits file photometry.fits
trig_time=`gethead ALEVT $file1 | sed -e 's/-//g' -e 's/://g'`
[ "$trig_time" ] || trig_time=`gethead ALALT $file1 | sed -e 's/-//g' -e 's/://g'`
[ "$trig_time" ] || trig_time=`gethead DATE-OBS stack_${cam}.fits | sed -e 's/-//g' -e 's/://g'`

tstart=
tstop=
if [ "$trig_time" ]; then
  gps0=`ut2gps.py $trig_time`
  ut1=`gethead DATE-OBS $file0 | sed -e 's/-//g' -e 's/://g'`
  ut2=`gethead DATE-OBE $file1 | sed -e 's/-//g' -e 's/://g'`
  tstart=`ut2gps.py $ut1 | awk '{printf("%.4f\n",($1-'$gps0')/3600)}'`
  tstop=`ut2gps.py $ut2 | awk '{printf("%.4f\n",($1-'$gps0')/3600)}'`
fi

cfilter=`gethead CFILTER stack_${cam}_dir/calibration.fits`
[ "$cfilter" ] || cfilter="GAIA-EDR3 G(AB)"

# dump the photometry to a master fits file
echo phot2fits_ddoti.py times_$file_list $ref $trig_time "$cfilter"
phot2fits_ddoti.py times_$file_list $ref $trig_time "$cfilter"

# make a data quality plot
echo zlf_plot_ddoti.py photometry.fits
zlf_plot_ddoti.py photometry.fits
mv photometry.fits.jpg zlf_plot.jpg

# make an astrometry plot
echo ast_err_plot.py
ast_err_plot.py

# make some lightcurves and ascii source lists
rm lc_*.* source_fitting.txt mp_phot.txt all_sources.txt new_phot.txt cat_only_phot.txt cal_phot.txt 2>/dev/null
echo ddoti_lc_plots.py photometry.fits $snr_cut
ddoti_lc_plots.py photometry.fits $snr_cut

for file in `ls new_phot.txt cat_only_phot.txt cal_phot.txt`; do
    head -1 $file > ${file}.tmp
    sed -n '2,$p' $file | sort -n -k 6 | sort -ns -k 12 >> ${file}.tmp
    mv ${file}.tmp $file
done

head -1 source_fitting.txt 2>/dev/null > source_fitting.txt.tmp
sed -n '2,$p' source_fitting.txt 2>/dev/null | sort -n -k 6 | sort -ns -k 17 >> source_fitting.txt.tmp
mv source_fitting.txt.tmp source_fitting.txt

c2file=../redux_colibri_AL/stack_AL.fits
if [ -f $c2file ]; then
    nx=`gethead NAXIS1 $c2file`
    ny=`gethead NAXIS2 $c2file`
    x0=`gethead CRPIX1 $c2file`
    y0=`gethead CRPIX2 $c2file`
    ra0=`gethead CRVAL1 $c2file`
    dec0=`gethead CRVAL2 $c2file`
    i00=`echo $nx $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
    j00=`echo $ny $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
    echo "$ra0 $dec0" > radec$$.txt
    radec2xy.py stack_${cam}.fits radec$$.txt > xy$$.txt
    x=`awk '{print $1;exit}' xy$$.txt`
    y=`awk '{print $2;exit}' xy$$.txt`
    i0=`echo $i00 $x $x0 | awk '{printf("%.0f\n",$1+$2-$3)}'`
    j0=`echo $j00 $y $y0 | awk '{printf("%.0f\n",$1+$2-$3)}'`
    rm radec$$.txt xy$$.txt 2>/dev/null
fi

# make the main jpg and postage stamps for the wegbpage
echo make_circles_ddoti.sh stack_${cam}.fits source_fitting.txt new_phot.txt cat_only_phot.txt all_sources.txt make_thumbs=yes thumb_size=$thumb_size thumb_edge=$thumb_edge thumb_scale=$thumb_scale logscale=$thumb_logscale invert_cmap=$thumb_invert make_fits_thumbs=$make_fits_thumbs i0=$i0 j0=$j0
make_circles_ddoti.sh stack_${cam}.fits source_fitting.txt new_phot.txt cat_only_phot.txt all_sources.txt make_thumbs=yes thumb_size=$thumb_size thumb_edge=$thumb_edge thumb_scale=$thumb_scale logscale=$thumb_logscale invert_cmap=$thumb_invert make_fits_thumbs=$make_fits_thumbs i0=$i0 j0=$j0

# start dumping stuff to web pages
# cal_phot_${cam}.html  :  indx>0
# new_phot_faint_${cam}.txt : indx<0 and snr<snr_cut
#
# some headers
#echo "#   id    RA           DEC        mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname" > hdr$$.txt
echo "#   id    RA           DEC                             mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname" > hdr$$.txt

#head -1 source_fitting.txt 2>/dev/null | sed -e 's/x,y//g' > hdr$$.txt
nfile=new_phot_${cam}.txt.summary.txt
vfile=variable_new_phot_${cam}.txt.summary.txt
sfile=spurious_new_phot_${cam}.txt.summary.txt
cp hdr$$.txt $nfile ; cp hdr$$.txt $vfile ; cp hdr$$.txt $sfile
rm hdr$$.txt

#
# this is how source id's are linked to thumbnail images
#
nx=`gethead NAXIS1 stack_${cam}.fits`
ny=`gethead NAXIS2 stack_${cam}.fits`
[ "$i0" ] || i0=`echo $nx $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
[ "$j0" ] || j0=`echo $ny $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
sethead TI0=$i0 TJ0=$j0 TS=$thumb_size TE=$thumb_edge PDIR=$pdir0 photometry.fits

fwhm0=`gethead FWHM stack_${cam}.fits`
sed -n '2,$p' source_fitting.txt 2>/dev/null | awk '{good=1; if($10>3*'$fwhm0') good=0; if ($13>10 && sqrt($11*$11)<0.1) good=0; if ($8<13) good=0; split($14,nepoch,"/"); if (nepoch[1]<0.75*nepoch[2] && $11>-$12) good=0; if ($15!=0 || $20>0) good=0; if (good==0) ofile="'$sfile'"; else if ($16==1) ofile="'$vfile'"; else ofile="'$nfile'"; x=int(($6-('$i0')-0.5)/'$thumb_size')*'$thumb_size'+'$i0';y=int(($7-('$j0')-0.5)/'$thumb_size')*'$thumb_size'+'$j0'; printf("<A HREF=\"'$pdir0'/thumb_%d_%d_.html\" title=\"Source %d, w=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" target=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %9s %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",x,y,$1,$8,$9,$20,$21,$19,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22)>>ofile}'

grep '#' new_phot.txt 2>/dev/null | tr 'x,y' ' ' > new_phot_faint_${cam}.txt
grep -v '#' new_phot.txt 2>/dev/null |  awk '{x=int(($6-('$i0')-0.5)/'$thumb_size')*'$thumb_size'+'$i0';y=int(($7-('$j0')-0.5)/'$thumb_size')*'$thumb_size'+'$j0'; printf("<A HREF=\"thumb_%d_%d_.html\" target=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %12.6f %12.6f %8.2f %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",x,y,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18)}' >> new_phot_faint_${cam}.txt

grep '#' mp_phot.txt 2>/dev/null | tr 'x,y' ' ' > mp_phot_${cam}.txt
grep -v '#' mp_phot.txt 2>/dev/null |  awk '{if($9<0) {x=int(($4-('$i0')-0.5)/'$thumb_size')*'$thumb_size'+'$i0';y=int(($5-('$j0')-0.5)/'$thumb_size')*'$thumb_size'+'$j0'; printf("<A HREF=\"thumb_%d_%d_.html\" target=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s %12.6f %12.6f %8.2f %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %s\n",x,y,$1,$2,$3,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16)}}' | sort -k 13 >> mp_phot_${cam}.txt

grep '#' cal_phot.txt 2>/dev/null | tr 'x,y' ' ' > cal_phot_${cam}.txt
grep -v '#' cal_phot.txt 2>/dev/null | awk '{x=int(($6-('$i0')-0.5)/'$thumb_size')*'$thumb_size'+'$i0';y=int(($7-('$j0')-0.5)/'$thumb_size')*'$thumb_size'+'$j0'; printf("<A HREF=\"thumb_%d_%d_.html\" target=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %12.6f %12.6f %8.2f %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",x,y,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18)}' >> cal_phot_${cam}.txt

grep '#' cat_only_phot.txt 2>/dev/null | tr 'x,y' ' ' > cat_only_phot_${cam}.txt
grep -v '#' cat_only_phot.txt 2>/dev/null | awk '{x=int(($6-('$i0')-0.5)/'$thumb_size')*'$thumb_size'+'$i0';y=int(($7-('$j0')-0.5)/'$thumb_size')*'$thumb_size'+'$j0'; printf("<A HREF=\"thumb_%d_%d_.html\" target=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %12.6f %12.6f %8.2f %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",x,y,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18)}' >> cat_only_phot_${cam}.txt

#
# compress the stacks for storage purposes
#
rm stack_${cam}.fits.fz stack_${cam}.wt.fits.fz 2>/dev/null
fpack stack_${cam}.fits &
fpack stack_${cam}.wt.fits &
wait

if [ -f stack_${cam}.diff.fits ]; then
    fpack stack_${cam}.diff.fits &
    fpack stack_${cam}.diff.wt.fits &
    wait
fi

pdir=${savedir}/$pdir0
[ -d "${pdir}" ] && rm -r $pdir
mkdir $pdir

cp ${REDUX_BASE_DIR}/jquery/jquery.maphilight.js $pdir
cp ${REDUX_BASE_DIR}/jquery/jQueryRotate.js $pdir

#
#make a webpage
#
exptime=`gethead EXPTIME stack_${cam}.fits | awk '{printf("%.0f\n",$1)}'`
fwhm=`gethead FWHM stack_${cam}.fits`
maglim=`gethead MAGLIM stack_${cam}.fits | awk '{printf("%.2f\n",$1)}'`
maglim5=`echo $maglim | awk '{printf("%.2f\n",$1+0.752575)}'`
maglim3=`echo $maglim | awk '{printf("%.2f\n",$1+1.307197)}'`
magzero=`gethead MAGZERO stack_${cam}.fits`

ps1file=`ls ps1_stack*.fits | egrep -v 'wt|thumb'`
maglimp=
[ -f "$ps1file" ] && maglimp=`gethead MAGLIM $ps1file | awk '{printf("%.2f\n",$1)}'`
if [ "$maglimp" ]; then
    maglimp5=`echo $maglimp | awk '{printf("%.2f\n",$1+0.752575)}'`
    maglimp3=`echo $maglimp | awk '{printf("%.2f\n",$1+1.307197)}'`
fi

detlim=
[ "$cam" = "AL" ] && detlim=`gethead DETLIM catalog_matches.fits`

maglims=
if [ -f subtraction_matches.fits ]; then
    maglims=`gethead MAGLIM subtraction_matches.fits`
    maglims5=`echo $maglims | awk '{printf("%.2f\n",$1+0.752575)}'`
    maglims3=`echo $maglims | awk '{printf("%.2f\n",$1+1.307197)}'`
fi

echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux</TITLE></HEAD>
<script>
function diffImage(src) {
       document.getElementById(\"mainimage\").src = src;
}
</script>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+2\" COLOR=\"#006600\"><A HREF=\"stack_${dte1}_${tag}_${nfiles}.html\" TARGET=\"_blank\">COLIBRI $cam $filter</A> RA $grb_RA , Dec $grb_DEC [$nfiles Frame(s), $dte1 - $dte2]</FONT><BR>
<FONT SIZE=\"+1\"> Size: $nx x $ny &nbsp; &nbsp; Exposure: ${exptime}s &nbsp; &nbsp; Time Since Trigger: $tstart - $tstop (hours) <BR>
  10-sigma limit: $maglim ($maglim5 5-sigma, $maglim3 3-sigma) <BR>
  Image Subtraction 10-sigma limit: $maglims ($maglims5 5-sigma, $maglims3 3-sigma) <BR>
  PS1 10-sigma limit: $maglimp ($maglimp5 5-sigma, $maglimp3 3-sigma) <BR>
  Zero-Point: $magzero &nbsp; &nbsp; Calibration Scatter: $calscatt &nbsp; &nbsp; FWHM: $fwhm pix </FONT><BR>" > index.html

cp index.html index_nocircles.html
cp index.html index_faint.html
cp index.html index_mp.html
cp index.html cal_phot_${cam}.html
cp index.html cat_only_phot_${cam}.html

echo "<IMG SRC=\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\" ID=\"mainimage\" USEMAP=\"#stack_${cam}.map\">" >> index.html
sed \$d stack_${cam}.map | sed -e "s/thumb_/${pdir0}\/thumb_/g" >> index.html

ps1file_jpeg=`ls ps1_stack_*.jpg`
diff_img=stack_${cam}.diff.jpg

echo "</map><BR>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\")>Standard view</Button> &nbsp; &nbsp;
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.jpg\")>Image without circles</BUTTON> &nbsp; &nbsp;" >> index.html
[ -f "$ps1file_jpeg" ] && echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/${ps1file_jpeg}\")>PS1 image</BUTTON> &nbsp; &nbsp;" >> index.html
[ -f "$diff_img" ] && echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.diff.jpg\")>${PS1}Difference image</BUTTON><P>" >> index.html

# see if there's source error region and report sources in there first 
if [ "$pos_err0" ]; then
    radec_sex=`radec2sex.py $ra0 $dec0`
    echo "<P><FONT SIZE=\"+1\" COLOR=\"#006600\"> Note: there is a $pos_err_as arcsec ($pos_err pixel) radius error region centered at: $radec_sex ($ra0 , $dec0) </FONT><BR><HR>" >> index.html
fi

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">$cfilter Band Photometry (Faint (SNR<$snr_cut) Uncatalogued Sources):</FONT><BR><P><PRE>" >> index_faint.html
cat new_phot_faint_${cam}.txt >> index_faint.html
echo "</PRE></BODY></HTML>" >> index_faint.html

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">$cfilter Band Photometry (Minor Planets):</FONT><BR><P><PRE>" >> index_mp.html
cat mp_phot_${cam}.txt >> index_mp.html
echo "</PRE></BODY></HTML>" >> index_mp.html

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">$cfilter Band Calibration Photometry (Bright (SNR>=$snr_cut) Catalogued Sources):</FONT><BR><P><PRE>" >> cal_phot_${cam}.html
cat cal_phot_${cam}.txt >> cal_phot_${cam}.html
echo "</PRE></BODY></HTML>" >> cal_phot_${cam}.html

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">$cfilter Band Photometry (Faint Catalogued Sources):</FONT><BR><P><PRE>" >> cat_only_phot_${cam}.html
cat cat_only_phot_${cam}.txt >> cat_only_phot_${cam}.html
echo "</PRE></BODY></HTML>" >> cat_only_phot_${cam}.html

echo "<P>Notes: high_pm=1 (USNO high proper motion star) high_pm=-1 (skybot minor planet) fading?=1 (fading source) fading?=-1 (brightening source) <BR><A HREF=\"${pdir0}/index_mp.html\" target=\"_blank\">minor planet matches</A><P>" >> index.html

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">$cfilter Band Photometry (Bright/Fading (SNR>=$snr_cut) Uncatalogued Sources):</FONT><BR><P><PRE>" >> index.html
cat $vfile >> index.html

echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Bright/Non-Fading (SNR>=$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> index.html
cat $nfile >> index.html

echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Bright/Questionable (SNR>=$snr_cut) Sources:</FONT><BR><P><PRE>" >> index.html
cat $sfile >> index.html

phot_diam=`gethead PHOTDIAM stack_${cam}.fits`
echo "</PRE><! sourcelists_done>
<A HREF=\"${pdir0}/index_faint.html\" target=\"_blank\">Faint New Source Photometry.</A> &nbsp;
<A HREF=\"${pdir0}/cal_phot_${dte1}_${tag}.html\" target=\"_blank\">Calibration Source Photometry.</A> &nbsp;
<A HREF=\"${pdir0}/cat_only_phot_${dte1}_${tag}.html\" target=\"_blank\">Faint Catalogued Source Photometry.</A>
<HR>
Using psf fitting diameter: $phot_diam (pixels)<P> Plots: " >> index.html
[ -f "zlf_plot.jpg" ] && echo "<A HREF=\"${pdir0}/zlf_plot.jpg\" TARGET=\"_blank\"> Sensitivity </A> &nbsp;" >> index.html
if [ -f psf_grid.fits ]; then
    fits2jpg.py psf_grid.fits 512 linvert
    [ -f "psf_grid.jpg" ] && echo "<A HREF=\"${pdir0}/psf_grid.jpg\" TARGET=\"_blank\"> PSF </A> &nbsp;" >> index.html

fi
[ -f "boxplot.jpg" ] && echo "<A HREF=\"${pdir0}/boxplot.jpg\" TARGET=\"_blank\"> Astrometry </A> &nbsp;" >> index.html
[ -f "scatter_plot.jpg" ] && echo "<A HREF=\"${pdir0}/scatter_plot.jpg\" TARGET=\"_blank\"> Lightcurve Scatter </A> &nbsp;" >> index.html
[ -f "stack_${cam}_dir/calibrated_catalog.jpg" ] && echo "<A HREF=\"${pdir0}/calibrated_catalog.jpg\" TARGET=\"_blank\"> Calibration Scatter </A> &nbsp;" >> index.html
[ -f "stack_${cam}_dir/calibrated_catalog_ra.jpg" ] && echo "<A HREF=\"${pdir0}/calibrated_catalog_ra.jpg\" TARGET=\"_blank\"> Calibration Scatter in RA </A> &nbsp;" >> index.html
[ -f "stack_${cam}_dir/calibrated_catalog_dec.jpg" ] && echo "<A HREF=\"${pdir0}/calibrated_catalog_dec.jpg\" TARGET=\"_blank\"> Calibration Scatter in DEC </A> &nbsp;" >> index.html
[ -f "distortion.jpg" ] && echo "<A HREF=\"${pdir0}/distortion.jpg\" TARGET=\"_blank\"> Residual Distortion </A> &nbsp;" >> index.html
[ -f "fwhm.jpg" ] && echo "<A HREF=\"${pdir0}/fwhm.jpg\" TARGET=\"_blank\"> Image FWHM </A> &nbsp;" >> index.html

cp $PYTHONPATH/photometry_ddoti.fits.txt $pdir
echo "<P><A HREF=\"${pdir0}/photometry_${dte1}_${tag}.fits\" TARGET=\"_blank\">All Source Photometry</A> (<A HREF=\"${pdir0}/photometry_ddoti.fits.txt\" TARGET=\"_blank\">Usage Description</A>).<BR>" >> index.html

dist_jpg=`ls *_dist.jpg 2>/dev/null | head -1`
[ -f "$dist_jpg" ] && echo "<A HREF=\"${pdir0}/stack_${dte1}_${tag}_dist.jpg\" target=\"_blank\">Distortion Correction</A>" >> index.html

task_tot=`date +%s | awk '{printf("Total Execution Time for All Tasks: %.2f minutes\n",($1-'$t00')/60.)}'`
upt=`uptime`
echo $upt

dte=`date -u`
echo ${savedir}/stack_${dte1}_${tag}_${nfiles}.log > logfile_name.txt
echo "
<HR WIDTH=\"100%\"> $task_tot
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.log\" target=\"_blank\">logfile</A>
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.fits.fz\" target=\"_blank\">stackfile</A>
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.wt.fits.fz\" target=\"_blank\">weightfile</A>" >> index.html
if [ -f stack_${cam}.diff.fits.fz ]; then
    echo "<A HREF=\"stack_${dte1}_${tag}_${nfiles}.diff.fits.fz\" target=\"_blank\">difffile</A>
          <A HREF=\"stack_${dte1}_${tag}_${nfiles}.diff.wt.fits.fz\" target=\"_blank\">diffweightfile</A>" >> index.html
fi
echo "<BR> $upt <BR> Last Updated: $dte (natbutler@asu.edu)</BODY></HTML>" >> index.html
echo "</BODY></HTML>" >> index_nocircles.html

#
# ship out all the summary files
#
cp index.html ${savedir}/stack_${dte1}_${tag}_${nfiles}.html
cp index_faint.html index_mp.html index_nocircles.html $pdir
rm ${savedir}/current_${cam}.html 2>/dev/null
ln -s ${savedir}/stack_${dte1}_${tag}_${nfiles}.html ${savedir}/current_${cam}.html

cp cal_phot_${cam}.html ${pdir}/cal_phot_${dte1}_${tag}.html
cp cat_only_phot_${cam}.html ${pdir}/cat_only_phot_${dte1}_${tag}.html
[ -f "$dist_jpg" ] && cp $dist_jpg ${pdir}/stack_${dte1}_${tag}_dist.jpg
cp new_phot_${cam}.txt.summary.txt ${pdir}/new_phot_${dte1}_${tag}_summary.txt
mv lc_*.jpg lc_*.txt lc_*.html $pdir 2>/dev/null
mv ps1_thumb_*_*_*_*_.jpg diff_thumb_*_*_*_*_.jpg thumb_*_*_*_*_.jpg thumb_*_*_*_*nocircles.jpg diff_thumb_*_*_.html thumb_*_*_.html thumb_*_*_nocircles.html $pdir 2>/dev/null
mv ps1_thumb_*_*_.jpg diff_thumb_*_*_.jpg thumb_*_*_.jpg $pdir 2>/dev/null
mv stack_${cam}_thumb.fits ${pdir}/stack_${dte1}_${tag}_thumb.fits
mv stack_${cam}_thumb.wt.fits ${pdir}/stack_${dte1}_${tag}_thumb.wt.fits
cp stack_${cam}.jpg ${pdir}/stack_${dte1}_${tag}.jpg
[ -f stack_${cam}.diff.jpg ] && mv stack_${cam}.diff.jpg ${pdir}/stack_${dte1}_${tag}.diff.jpg
[ -f $ps1file_jpeg ] && cp $ps1file_jpeg ${pdir}/$ps1file_jpeg
mv stack_${cam}_circles.jpg ${pdir}/stack_${dte1}_${tag}_circles.jpg
mv stack_${cam}.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.fits.fz
mv stack_${cam}.wt.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.wt.fits.fz
if [ -f stack_${cam}.diff.fits.fz ]; then
    mv stack_${cam}.diff.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.diff.fits.fz
    mv stack_${cam}.diff.wt.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.diff.wt.fits.fz
fi
cp scatter_plot.jpg zlf_plot.jpg psf_grid.jpg boxplot.jpg distortion.jpg fwhm.jpg stack_${cam}_dir/calibrated_catalog*.jpg $pdir 2>/dev/null
cp photometry.fits ${pdir}/photometry_${dte1}_${tag}.fits

# cleanup
rm thumb.batch phot$$.tmp0 phot$$.tmp1 thumb_*_*_*_*_.map stack_*.map usno_xy.txt zlf_list$$.txt cal_phot_${cam}.tmp 2>/dev/null
