#!/bin/bash

file_list=$1
shift

savedir=`pwd`
make_frame_jpegs=no
logfile=
catdir=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

file0=`sort $file_list | head -1`
cam=`gethead CCD_NAME $file0 | cut -c1-2`
[ "$cam" ] || cam=`basename $file0 | cut -c17-18`
[ "$filter" ] || filter=`gethead FILTER stack_${cam}.fits`
[ "$filter" ] || filter=r

# cfilter is the calibration filter
cfilter="${filter}-PS1-DDRAGO"

dte1=`echo $file0 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`
file1=`sort $file_list | tail -1`
dte2=`echo $file1 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`

tag=${cam}_$filter
grb_RA=`gethead CRVAL1 stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
grb_DEC=`gethead CRVAL2 stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`

# make a jpg and one with circles
[ -f stack_${cam}.jpg ] && rm stack_${cam}.jpg
file=stack_${cam}.fits
x0=`gethead CRPIX1 $file`; y0=`gethead CRPIX2 $file`
dx0=`gethead NAXIS1 $file`; dy0=`gethead NAXIS2 $file`
radec2xy.py stack_${cam}.fits usno_radec.fits | awk '{if($1>=1 && $1<='$dx0' && $2>=1 && $2<='$dy0') printf("%f %f %.0f_%.1f\n",$1,$2,$3,$4)}' > usno_xy.txt
echo make_circles.sh stack_${cam}.fits stack_${cam}_radec.txt calfile=stack_${cam}_radec0.txt
make_circles.sh stack_${cam}.fits stack_${cam}_radec.txt calfile=stack_${cam}_radec0.txt

nfiles=`cat $file_list | wc -l`
pdir0=photometry_${dte1}_${tag}_${nfiles}
pdir=${savedir}/$pdir0
if [ -d $pdir ]; then
    rm -r ${pdir}/*
else
    mkdir $pdir
fi

if [ "$make_frame_jpegs" = "yes" ]; then
    for file in `cat $file_list`; do
        base=${file%'.fits'}
        echo fits2jpg.py $file 1024 linvert
        fits2jpg.py $file 1024 linvert &
        gonogo
    done
    wait; go_iter=0
fi

# get a ps1 or dss frame
ps1file=ps1_stack_${filter}.fits
diff_img=stack_${cam}.diff.fits
diff_wimg=stack_${cam}.diff.wt.fits
rm $diff_img $diff_wimg 2>/dev/null
echo get_ps1.sh stack_${cam}.fits outfile=$ps1file catdir=$catdir
get_ps1.sh stack_${cam}.fits outfile=$ps1file catdir=$catdir
ps1file_jpeg=ps1_stack_${filter}.jpg
[ -f "$ps1file" ] || get_dss.sh stack_${cam}.fits dss_base_file=${catdir}/dss0.fits

exptime=`gethead EXPTIME stack_${cam}.fits | awk '{printf("%.1f\n",$1)}'`
fwhm=`gethead FWHM stack_${cam}.fits`
maglim=`gethead MAGLIM stack_${cam}.fits`
magzero=`gethead MAGZERO stack_${cam}.fits`
pos_err0=`gethead ALUN stack_${cam}.fits`

#make a difference image
if [ -f "$diff_img" ]; then
    echo "Full frame difference image created."
elif [ "$pos_err0" ]; then
    rm stack_${cam}.diff.fits 2>/dev/null
    sz=`echo $pos_err0 | awk '{printf("%.0f\n", 2*$1*3600/0.25)}'`
    [ "$sz" -gt 2048 ] && sz=2048
    [ "$sz" -lt 1024 ] && sz=1024
    echo colibri_imsubPS1.sh stack_${cam}.fits sz=$sz fwhm=$fwhm catdir=$catdir
    colibri_imsubPS1.sh stack_${cam}.fits sz=$sz fwhm=$fwhm catdir=$catdir
    [ -f "stack_${cam}.diff.fits" ] && fits2jpg.py stack_${cam}.diff.fits 1024 linvert
fi
[ -f "$diff_img" ] && PS1="PS1 "
if [ "$PS1" ]; then
    echo "Using PS1 for difference image analysis."
elif [ -f stack_${cam}a.fits -a -f stack_${cam}b.fits ]; then
    echo "Using DSS for difference image analysis."
    diff_img=stack_${cam}a.diff.fits
    rm $diff_img 2>/dev/null
    echo coatli_imsubAB.sh stack_${cam}a.fits stack_${cam}b.fits fwhm=$fwhm
    coatli_imsubAB.sh stack_${cam}a.fits stack_${cam}b.fits fwhm=$fwhm
    [ -f stack_${cam}a.diff.fits ] && fits2jpg.py stack_${cam}a.diff.fits 1024 linvert
fi

fits2jpg.py stack_${cam}a.fits 1024 linvert
fits2jpg.py stack_${cam}b.fits 1024 linvert
fits2jpg.py stack_${cam}.wt.fits 1024
maskfile=`cat maskfile_used.txt`
fits2jpg.py $maskfile 1024

nx=`gethead NAXIS1 stack_${cam}.fits`
ny=`gethead NAXIS2 stack_${cam}.fits`

#
#make a webpage
#
echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux</TITLE></HEAD>
<script>
function diffImage(src) {
       document.getElementById(\"mainimage\").src = src;
}
</script>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+2\" COLOR=\"#006600\">COLIBRI $cam $filter : RA $grb_RA , Dec $grb_DEC
[N=$nfiles Frame(s), $dte1 - $dte2]</FONT><BR>
<FONT SIZE=\"+1\">
Frame Size: $nx x $ny &nbsp; &nbsp; &nbsp; Exposure Time: $exptime seconds &nbsp; &nbsp; &nbsp;
10-sigma limiting mag: $maglim &nbsp; &nbsp; &nbsp;
Zero-Point: $magzero &nbsp; &nbsp; &nbsp;
FWHM: $fwhm pixels </FONT><P>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\")>Standard view</Button>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.jpg\")>Image without circles</BUTTON>" >> index.html
if [ -f "$ps1file_jpeg" ]; then
    echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/${ps1file_jpeg}\")>PS1 image</BUTTON>" >> index.html
elif [ -f dss.jpg ]; then
    echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/dss.jpg\")>DSS image</BUTTON>" >> index.html
fi
echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.wt.jpg\")>Weight image</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}_mask.jpg\")>Mask image</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}a.jpg\")>First N/2 frames</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}b.jpg\")>Second N/2 frames</BUTTON>" >> index.html
[ -f "$diff_img" ] && echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.diff.jpg\")>${PS1}Difference image</BUTTON>" >> index.html

echo "<BR><IMG SRC=\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\" ID=\"mainimage\" USEMAP=\"#stack_${cam}_circles.map\">" >> index.html
sed \$d stack_${cam}_circles.map | sed -e "s/lc_/${pdir0}\/lc_/g" >> index.html

echo "</map><BR>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\")>Standard view</Button>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.jpg\")>Image without circles</BUTTON>" >> index.html
if [ -f "$ps1file_jpeg" ]; then
    echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/${ps1file_jpeg}\")>PS1 image</BUTTON>" >> index.html
elif [ -f dss.jpg ]; then
    echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/dss.jpg\")>DSS image</BUTTON>" >> index.html
fi
echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.wt.jpg\")>Weight image</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}_mask.jpg\")>Mask image</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}a.jpg\")>First N/2 frames</BUTTON>
<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}b.jpg\")>Second N/2 frames</BUTTON>" >> index.html
[ -f "$diff_img" ] && echo "<BUTTON ONCLICK=diffImage(\"${pdir0}/stack_${dte1}_${tag}.diff.jpg\")>${PS1}Difference image</BUTTON><P>" >> index.html


echo "<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux</TITLE></HEAD><BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+1\" COLOR=\"#006600\">Sources not in USNO-B1(<A HREF=\"source_fitting_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> Variability plot </A>):</FONT><BR><P><PRE>
#ID        RA          Dec           X            Y           mag         dmag            FWHM    Dist.  Flags" > new_sources.html

#
# go from radec.txt and radec0.txt to radec.tmp to radec0.tmp
#
# start by adding distance column

grep '#' stack_${cam}_radec.txt > stack_${cam}_radec.tmp
grep '#' stack_${cam}_radec0.txt > stack_${cam}_radec0.tmp
if [ "$pos_err0" ]; then
    pos_err1=`echo $pos_err0 | awk '{s=3600.*$1; printf("%.1f\n",sqrt(1+s*s))}'`
    sx0=`gethead SEXX0 stack_${cam}.fits`
    sy0=`gethead SEXY0 stack_${cam}.fits`
    ps=`gethead CD1_1 CD1_2 stack_${cam}.fits | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600.)}'`

    grep -v '#' stack_${cam}_radec.txt | awk '{x=$8;y=$9; dx=x-'$sx0'; dy=y-'$sy0'; dis=sqrt(dx*dx+dy*dy)*'$ps'; printf("%s %8.1f\n",$0,dis)}' >> stack_${cam}_radec.tmp
    grep -v '#' stack_${cam}_radec0.txt | awk '{x=$8;y=$9; dx=x-'$sx0'; dy=y-'$sy0'; dis=sqrt(dx*dx+dy*dy)*'$ps'; printf("%s %8.1f\n",$0,dis)}' >> stack_${cam}_radec0.tmp
else
    grep -v '#' stack_${cam}_radec.txt | awk '{printf("%s %8.1f\n",$0,0.0)}' >> stack_${cam}_radec.tmp
    grep -v '#' stack_${cam}_radec0.txt | awk '{printf("%s %8.1f\n",$0,0.0)}' >> stack_${cam}_radec0.tmp
fi

# flag fading sources among the uncatalogued sources
if [ -s fading_sources.txt ]; then
    awk 'BEGIN{print ""; n=-1}{for(i=n;i>$1;i--) print ""; print "fading_source"; n=$1-1}' fading_sources.txt > faders$$.txt
    paste -d " " stack_${cam}_radec.tmp faders$$.txt > stack_${cam}_radec.tmp1
    mv stack_${cam}_radec.tmp1 stack_${cam}_radec.tmp
    rm faders$$.txt
fi

# flag ps1_source_not_in_usno
#matchfile=stack_${cam}_dir/catalog.fits.ps1_matched.txt
#[ -s "$matchfile" ] || matchfile=stack_${cam}_radec.txt.matched.txt
matchfile=stack_${cam}_radec.txt.matched.txt
if [ -s "$matchfile" ]; then
    awk '{if($1==1) print "ps1_source_not_in_usno"; else print ""}' $matchfile > matched$$.txt
    paste -d " " stack_${cam}_radec.tmp matched$$.txt > stack_${cam}_radec.tmp1
    mv stack_${cam}_radec.tmp1 stack_${cam}_radec.tmp
    rm matched$$.txt
fi

# flag sources detected in image subtraction
isub_file=stack_${cam}.diff_dir/stack_${cam}.diff_imsub.txt
if [ -s "$isub_file" ]; then
    awk 'BEGIN{print ""; n=-1}{for(i=n;i>$1;i--) print ""; print "subtraction_detected"; n=$1-1}' $isub_file > imsubers$$.txt
    paste -d " " stack_${cam}_radec.tmp imsubers$$.txt > stack_${cam}_radec.tmp1
    mv stack_${cam}_radec.tmp1 stack_${cam}_radec.tmp
    rm imsubers$$.txt
fi

#
# need to fix these slow bits
#

nnew=`grep -v '#' stack_${cam}_radec.tmp | wc -l`
ls lc_dir/lc_*txt | sed -e 's/lc_dir\/lc_//g' -e 's/\.txt//g' | sort -nr |\
       awk 'BEGIN{n=-1}{for(i=n;i>$1;i--) print i; printf("<A|HREF=\"lc_%d.jpg\"|TARGET=\"_blank\">%d</A><A|HREF=\"lc_%d.txt\"|TARGET=\"_blank\">t</A>\n",$1,$1,$1); n=$1-1}END{for(i=n;i>=-'$nnew';i--) print i}' > tag_radec.txt

#grep -v '#' stack_${cam}_radec.tmp | sort -n -k 16 | awk '{if (system("[ -f lc_dir/lc_"$15".jpg ]")) printf("%d\n",$15); else printf("<A HREF=\"lc_%d.jpg\" TARGET=\"_blank\">%d</A><A HREF=\"lc_%d.txt\" TARGET=\"_blank\">t</A>\n",$15,$15,$15)}' > tag_radec.txt
#grep -v '#' stack_${cam}_radec.tmp | sort -n -k 16 | awk '{s=""; if(NF>16) {for(i=17;i<=NF;i++) s=s" "$i}; printf("%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f %8.1f %s\n",$1,$2,$8,$9,$3,$4,$7,$16,s)}' > rest_radec.txt
grep -v '#' stack_${cam}_radec.tmp | awk '{s=""; if(NF>16) {for(i=17;i<=NF;i++) s=s" "$i}; printf("%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f %8.1f %s\n",$1,$2,$8,$9,$3,$4,$7,$16,s)}' > rest_radec.txt

if [ "$pos_err0" ]; then
    paste tag_radec.txt rest_radec.txt | sort -n -k 9 | sed -e 's/|/ /g' >> new_sources.html
else
    paste tag_radec.txt rest_radec.txt | sed -e 's/|/ /g' >> new_sources.html
fi
echo "</PRE></BODY></HTML>" >> new_sources.html

rm stack_${cam}_radec.tmp

echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux (raw)</TITLE></HEAD>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+1\" COLOR=\"#006600\">Sources not in USNO-B1 (uncalibrated photometry):</FONT><BR><P><PRE>
#ID        RA          Dec           X            Y           mag         dmag            FWHM" > index_raw.html
grep -v '#' stack_${cam}_radec_raw.txt | awk '{printf("%-10s %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f\n",$15,$1,$2,$8,$9,$3,$4,$7)}' >> index_raw.html
echo "</PRE></BODY></HTML>" >> index_raw.html

echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux (raw)</TITLE></HEAD>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+1\" COLOR=\"#006600\">$cfilter Band Photometry (Catalogued Sources,
<A HREF=\"calplot_${dte1}_${tag}a.jpg\" TARGET=\"_blank\"> Calibration plot </A>,
<A HREF=\"scatter_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> Photometric Scatter plot </A>,
<A HREF=\"psf_grid_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> PSF Grid plot </A>,
<A HREF=\"zlf_plot.jpg\" TARGET=\"_blank\"> Sensitivity plot </A>):</FONT><BR><P><PRE>
#ID        RA          Dec           X            Y           mag         dmag            FWHM    Dist.   Flags" > catalogued_sources.html

grep '#' stack_${cam}_radec0.tmp > stack_${cam}_radec0.tmp1
grep -v '#' stack_${cam}_radec0.tmp | awk '{if($15<0) print $0,"usno_high_pm"; else print}' >> stack_${cam}_radec0.tmp1
mv stack_${cam}_radec0.tmp1 stack_${cam}_radec0.tmp

# flag usno_source_not_in_ps1
if [ -s stack_${cam}_radec0.txt.matched.txt ]; then
    awk '{if($1==0) print "usno_source_not_in_ps1"; else print ""}' stack_${cam}_radec0.txt.matched.txt > notmatched$$.txt
    paste -d " " stack_${cam}_radec0.tmp notmatched$$.txt > stack_${cam}_radec0.tmp1
    mv stack_${cam}_radec0.tmp1 stack_${cam}_radec0.tmp
    rm notmatched$$.txt
fi

#
# don't need links for catalogued sources
#

nold=`grep -v '#' stack_${cam}_radec0.tmp | wc -l`
grep -v '#' stack_${cam}_radec0.tmp | sort -n -k 16 | awk '{if (system("[ -f lc_"$15".jpg ]")) printf("%d\n",$15); else printf("<A HREF=\"lc_%d.jpg\" TARGET=\"_blank\">%d</A><A HREF=\"lc_%d.txt\" TARGET=\"_blank\">t</A>\n",$15,$15,$15)}' > tag_radec0.txt
grep -v '#' stack_${cam}_radec0.tmp | sort -n -k 16 | awk '{s=""; if(NF>16) {for(i=17;i<=NF;i++) s=s" "$i}; printf("%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f %8.1f %s\n",$1,$2,$8,$9,$3,$4,$7,$16,s)}' > rest_radec0.txt
paste tag_radec0.txt rest_radec0.txt >> catalogued_sources.html
echo "</PRE></BODY></HTML>" >> catalogued_sources.html
rm stack_${cam}_radec0.tmp

# if there's source error region and report sources in there first
if [ "$pos_err0" ]; then
    cat new_sources.html catalogued_sources.html | grep -v '#' > index.tmp1
    pos_err=`echo $pos_err0 | awk '{printf("%.1f\n",3600.*$1)}'`
    pos_err1=`echo $pos_err0 | awk '{s=3600.*$1; printf("%.1f\n",sqrt(1+s*s))}'`
    ra0=`gethead ALRA stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead ALDE stack_${cam}.fits | awk '{printf("%.6f\n",$1)}'`
    echo "<FONT SIZE=\"+2\" COLOR=\"#006600\"> Sources within $pos_err arcsec from center: $ra0 , $dec0 </FONT><BR><P><PRE>" >> index.html
    grep '#ID' catalogued_sources.html | awk '{print;exit}' >> index.html
    grep -v '#' index.tmp1 | awk '{if($2~/HREF/) dis=$13; else dis=$9; if (dis<'$pos_err1' && NF>8) print}' | sed -e "s/lc_/${pdir0}\/lc_/g" >> index.html
    echo "</PRE><P><HR>" >> index.html
    rm index.tmp1
fi

echo "<FONT SIZE=\"+1\" COLOR=\"#006600\"><A HREF=\"${pdir0}/new_sources.html\" TARGET=\"_blank\">$nnew Sources not in USNO-B1</A>
(<A HREF=\"${pdir0}/source_fitting_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> Variability plot </A>)</FONT><BR><P>
<FONT SIZE=\"+1\" COLOR=\"#006600\"><A HREF=\"${pdir0}/catalogued_sources.html\" TARGET=\"_blank\">$nold Catalogued Sources</A>
($cfilter Band Photometry, <A HREF=\"${pdir0}/calplot_${dte1}_${tag}a.jpg\" TARGET=\"_blank\"> Calibration plot </A>,
<A HREF=\"${pdir0}/scatter_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> Photometric Scatter plot </A>,
<A HREF=\"${pdir0}/psf_grid_${dte1}_${tag}.jpg\" TARGET=\"_blank\"> PSF Grid plot </A>,
<A HREF=\"${pdir0}/zlf_plot.jpg\" TARGET=\"_blank\"> Sensitivity plot </A>)</FONT><BR><P>" >> index.html

echo "</PRE><HR><FONT SIZE=\"+1\" COLOR=\"#006600\">$cfilter Band Photometry (Catalogued Sources, uncalibrated photometry):</FONT><BR><P><PRE>
#ID        RA          Dec           X            Y           mag         dmag            FWHM" >> index_raw.html
#grep -v '#' stack_${cam}_radec0_raw.txt | awk '{if (system("[ -f lc_dir/lc_"$NF".jpg ]")) printf("%d\n",$NF); else printf("<A HREF=\"lc_%d.jpg\" TARGET=\"_blank\">%d</A><A HREF=\"lc_%d.txt\" TARGET=\"_blank\">t</A>\n",$NF,$NF,$NF)}' > tag_radec0.txt
#grep -v '#' stack_${cam}_radec0_raw.txt | awk '{printf("%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f\n",$1,$2,$8,$9,$3,$4,$7)}' > rest_radec0.txt
grep -v '#' stack_${cam}_radec0_raw.txt | awk '{printf("%-10s %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %8.2f\n",$15,$1,$2,$8,$9,$3,$4,$7)}' >> index_raw.html
#paste tag_radec0.txt rest_radec0.txt >> index_raw.html
echo "</PRE></BODY></HTML>" >> index_raw.html

cp $PYTHONPATH/photometry.fits.txt $pdir
dmag0=`gethead CALERR stack_${cam}.fits`
echo "<A HREF=\"${pdir0}/photometry_${dte1}_${tag}.fits\" TARGET=\"_blank\">All Source Photometry</A> (<A HREF=\"${pdir0}/photometry.fits.txt\" TARGET=\"_blank\">Usage Description</A>).<BR>
Note: all photometric errors include a $dmag0 calibration error term added in quadrature.
The uncalibrated photometry is <A HREF=\"${pdir0}/stack_${dte1}_${tag}_${nfiles}_raw.html\" TARGET=\"_blank\">here.</A>" >> index.html

if [ "$make_frame_jpegs" = "yes" ]; then
    echo "<HR> Individual Frame JPEGs:<BR>" >> index.html
    ls f20*.jpg | awk '{printf("<A HREF=\"'$pdir0'/%s\" TARGET=\"_blank\">%03d</A> ",$1,NR); if(NR%30==0) printf("<BR>\n")}END{printf("\n")}' >> index.html
fi

upt=`uptime`
echo $upt

dte=`date -u`
echo ${savedir}/stack_${dte1}_${tag}_${nfiles}.log > logfile_name.txt
echo "<HR WIDTH=\"100%\"> $task_tot 
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.log\" target=\"_blank\">logfile</A>
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.fits.fz\" target=\"_blank\">stackfile</A>
<A HREF=\"stack_${dte1}_${tag}_${nfiles}.wt.fits.fz\" target=\"_blank\">weightile</A><BR>
$upt <BR> Last Updated: $dte (natbutler@asu.edu)</BODY></HTML>" >> index.html

#
# ship out all the summary files
#
mv index.html ${savedir}/stack_${dte1}_${tag}_${nfiles}.html
[ -f index_raw.html ] && mv index_raw.html ${pdir}/stack_${dte1}_${tag}_${nfiles}_raw.html
[ -f new_sources.html ] && mv new_sources.html ${pdir}/new_sources.html
[ -f catalogued_sources.html ] && mv catalogued_sources.html ${pdir}/catalogued_sources.html
rm ${savedir}/current_${cam}.html 2>/dev/null
ln -s ${savedir}/stack_${dte1}_${tag}_${nfiles}.html ${savedir}/current_${cam}.html

cp stack_${cam}_radec.txt ${pdir}/stack_${dte1}_${tag}_radec.txt
cp stack_${cam}_radec0.txt ${pdir}/stack_${dte1}_${tag}_radec0.txt
cp lc_dir/lc_*.jpg lc_dir/lc_*.txt $pdir 2>/dev/null
mv f20*.jpg $pdir 2>/dev/null
mv stack_${cam}.jpg ${pdir}/stack_${dte1}_${tag}.jpg
[ -f stack_${cam}a.diff.jpg ] && mv stack_${cam}a.diff.jpg ${pdir}/stack_${dte1}_${tag}.diff.jpg
[ -f stack_${cam}.diff.jpg ] && mv stack_${cam}.diff.jpg ${pdir}/stack_${dte1}_${tag}.diff.jpg
mv stack_${cam}.wt.jpg ${pdir}/stack_${dte1}_${tag}.wt.jpg
mv stack_${cam}a.jpg ${pdir}/stack_${dte1}_${tag}a.jpg
mv stack_${cam}b.jpg ${pdir}/stack_${dte1}_${tag}b.jpg
mv ${maskfile%'.fits'}.jpg ${pdir}/stack_${dte1}_${tag}_mask.jpg
mv stack_${cam}_circles.jpg ${pdir}/stack_${dte1}_${tag}_circles.jpg
cp calplot.jpg ${pdir}/calplot_${dte1}_${tag}a.jpg
if [ -f psf_grid.fits ]; then
    fits2jpg.py psf_grid.fits 512 linvert
    cp psf_grid.jpg ${pdir}/psf_grid_${dte1}_${tag}.jpg
fi
cp scatter_plot.jpg ${pdir}/scatter_${dte1}_${tag}.jpg
cp source_fitting.txt.jpg ${pdir}/source_fitting_${dte1}_${tag}.jpg
cp zlf_plot.jpg $pdir
cp photometry.fits ${pdir}/photometry_${dte1}_${tag}.fits
[ -f dss.jpg ] && cp dss.jpg $pdir
[ -f "$ps1file_jpeg" ] && cp $ps1file_jpeg $pdir

rm stack_${cam}.fits.fz stack_${cam}.wt.fits.fz 2>/dev/null
fpack stack_${cam}.fits &
fpack stack_${cam}.wt.fits &
wait

mv stack_${cam}.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.fits.fz
mv stack_${cam}.wt.fits.fz ${savedir}/stack_${dte1}_${tag}_${nfiles}.wt.fits.fz

rm phot$$.tmp0 phot$$.tmp1 thumb_*.fits thumb_*.map stack_*.map usno_xy.txt 2>/dev/null
