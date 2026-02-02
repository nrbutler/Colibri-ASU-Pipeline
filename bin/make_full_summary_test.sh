#!/bin/bash

snr_cut=10
question=no
thumb_size=401
thumb_edge=30
thumb_scale=1

sz=1024
savedir=`pwd`

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

hdr0="#    id    RA           DEC                             mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname"
hdr="#    id   RA           DEC                                 mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname"

#
# collect sources from the individual webpages
#

rm bright_fading.txt bright_nonfading.txt bright_questionable.txt 2>/dev/null
nfiles0=0
filters=
for hfile in `ls redux_colibri_C?_*/index.html`; do
    photdir=`grep index_faint $hfile | awk '{print $2}' | sed -e 's/"/ /g' -e 's/\// /g' | awk '{print $2}'`
    filter=`echo $photdir | awk -F_ '{print $(NF-1)}'`
    nfiles=`echo $photdir | awk -F_ '{print $NF}'`
    nfiles0=`expr $nfiles0 + $nfiles`
    filters="${filters}$filter"
    awk '/Bright\/Fading/,/Bright\/Non-Fading/{print}' $hfile | egrep -v 'dslope|Fading' | sed -e "s/<A/${filter}<A/" >> bright_fading.txt
    awk '/Bright\/Non-Fading/,/Questionable/{print}' $hfile | egrep -v 'dslope|Bright' | sed -e "s/<A/${filter}<A/" >> bright_nonfading.txt
    [ "$question" = "yes" ] && awk '/Questionable/,/sourcelists_done/{print}' $hfile | egrep -v 'dslope|Bright|sourcelists' | sed -e "s/<A/${filter}<A/" >> bright_questionable.txt
    dir=`dirname $hfile`
    grep HREF ${dir}/index_faint.html | grep thumb | sed -e "s/thumb/${photdir}\/thumb/g" -e "s/<A/${filter}<A/" >> faint_sources.txt
done
sort -nr -k 5 faint_sources.txt > faint_sources.tmp
mv faint_sources.tmp faint_sources.txt

tag="${filters}_${nfiles0}"
pdir0=summary_$tag
rm -r $pdir0 2>/dev/null
rm new_sources 2>/dev/null
mkdir $pdir0
ln -s $pdir0 new_sources

# webpage to make
index=summary_${tag}.html

sz0=`echo $thumb_size $thumb_edge $thumb_scale | awk '{print ($1+2*$2)*$3}'`

here=`pwd`
dir=`basename $here | awk -F_ '{print $6"/"$7"/"$8"/"$9}'`
# make a webpage
echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>Redux Summary</TITLE></HEAD>
<script>
function diffImage(src) {
       document.getElementById(\"mainimage\").src = src;
}
</script>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+2\" COLOR=\"#006600\">Summary Chi2 Image for ${dir} &nbsp; - &nbsp; $nfiles0 files &nbsp; - &nbsp; $filters filters</FONT><BR>
<IMG SRC=\"$pdir0/stack_AL_circles.jpg\" ID=\"mainimage\" USEMAP=\"#summary.map\"><BR>
<BUTTON ONCLICK=diffImage(\"$pdir0/stack_AL_circles.jpg\")>Standard view</Button> &nbsp; &nbsp;
<BUTTON ONCLICK=diffImage(\"$pdir0/stack_AL.jpg\")>Image without circles</BUTTON> &nbsp; &nbsp;" > $index
ps1file_jpeg=ps1_stack.jpg
[ -f "redux_colibri_AL/$ps1file_jpeg" ] && echo "<BUTTON ONCLICK=diffImage(\"$pdir0/${ps1file_jpeg}\")>PS1 image</BUTTON> &nbsp; &nbsp;" >> $index

echo "<P>" >> $index
awk '/COLIBRI C/,/FWHM:/{print}' redux_colibri_C?_*/index.html >> $index

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">Bright/Fading (SNR>=$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> $index
echo "$hdr0" >> $index
awk '{id=$4; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_fading.txt | sed -e 's/,\.html/\.html/g' >> $index
echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Bright/Non-Fading (SNR>=$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> $index
echo "$hdr0" >> $index
awk '{id=$4; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_nonfading.txt | sed -e 's/,\.html/\.html/g' >> $index
if [ -s bright_questionable.txt ]; then
    echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Bright/Questionable (SNR>=$snr_cut) Sources:</FONT><BR><P><PRE>" >> $index
    echo "$hdr0" >> $index
    awk '{id=$4; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_questionable.txt | sed -e 's/,\.html/\.html/g' >> $index
fi
echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Faint (SNR<$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> $index
echo "$hdr" >> $index
awk '{id=$5; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' faint_sources.txt | sed -e 's/<\/A>\.html/\.html/g' >> $index

rm faint_sources.txt bright_fading.txt bright_nonfading.txt bright_questionable.txt 2>/dev/null

# retrieve the candidate id numbers
rm candidate_ids.txt 2>/dev/null
grep HREF $index | grep thumb | grep title | awk '{print $4}' | sed -e 's/,//g' > candidate_ids.txt
grep HREF $index | grep thumb | grep -v title | awk '{print $5}' | sed -e 's/<\/A>//' >> candidate_ids.txt
sort -nr -u candidate_ids.txt > ids$$.tmp
mv ids$$.tmp candidate_ids.txt

# compile source info
source_info.py candidate_ids.txt $snr_cut

fitsfile=redux_colibri_AL/stack_AL.fits
pos_err0=`gethead ALUN $fitsfile | awk '{printf("%.6f\n",$1)}'`
if [ "$pos_err0" ]; then
    ps=`gethead CD1_1 CD1_2 $fitsfile | awk '{printf("%.2f\n",sqrt($1*$1+$2*$2)*3600.)}'`
    ra0=`gethead -c ALRA $fitsfile | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead -c ALDE $fitsfile | awk '{printf("%.6f\n",$1)}'`
    pos_err=`echo $pos_err0 $ps | awk '{printf("%.1f\n",3600.*$1/$2)}'`
    echo $ra0 $dec0 > radec$$.txt
    radec2xy.py $fitsfile radec$$.txt noverify > xy$$.txt
    xpos=`awk '{print $1;exit}' xy$$.txt`
    ypos=`awk '{print $2;exit}' xy$$.txt`
    echo $xpos $ypos $pos_err > radec0_$$.txt
    rm xy$$.txt radec$$.txt
fi

nx=`gethead NAXIS1 $fitsfile`
ny=`gethead NAXIS2 $fitsfile`
sx=$sz
sy=`echo "$sz $nx $ny" | awk '{printf("%.0f\n",$1*$3/$2)}'`
facx=`echo $sx $nx | awk '{printf("%.6f\n",$2/$1)}'`
facy=`echo $sy $ny | awk '{printf("%.6f\n",$2/$1)}'`

greencircles=
if [ "$pos_err0" ]; then
    greencircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke green -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -draw '\''line 0,%.0f %.0f,%.0f'\'' -draw '\''line %.0f,0 %.0f,%.0f'\'' ",x,y,x,y+$3/'$facy',y,('$nx'-1)/'$facx',y,x,x,('$ny'-1)/'$facy')}' radec0_$$.txt`
fi

radeclist=radec$$.txt
grep Source $pdir0/source_*.txt | awk '{print $5,$6,$2,$2}' | tr '(,:)' ' ' > $radeclist

redcircles=
if [ -s "$radeclist" ]; then

    radec2xy.py $fitsfile $radeclist noverify > xy$$.txt
    redcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';n=$4; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -stroke black -fill green -draw '\''translate %.0f,%.0f text -5,-11 \"%.0f\"'\''",x,y,x,y+5,x,y,n)}' xy$$.txt`

    echo "<map name=\"summary.map\">" > summary.map
    awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$4; printf("<area shape =\"circle\" coords=\"%.0f,%.0f,5\" href=\"'${pdir0}'/source_%.0f.html\" target=\"_blank\" title=\"Source %.0f\">\n",x,y,id,id)}' xy$$.txt >> summary.map
    echo "</map>" >> summary.map

fi
rm $radeclist xy$$.txt radec0_$$.txt 2>/dev/null

if [ -f summary.map ]; then
    cat summary.map >> $index
    rm summary.map
fi
dte=`date -u`
echo "</PRE><P><HR WIDTH=\"100%\">Last Updated: $dte (natbutler@asu.edu)</BODY></HTML>" >> $index

cmd="convert -quality 75 redux_colibri_AL/stack_AL.jpg -pointsize 16 $redcircles $greencircles redux_colibri_AL/stack_AL_circles.jpg"
eval "$cmd"

for nid in `cat candidate_ids.txt`; do
    file=new_sources/source_${nid}.txt
    if [ -s $file ]; then
        ofile=new_sources/source_${nid}.html
        echo "<HTML>" > $ofile
        cat ${REDUX_BASE_DIR}/jquery/hdr.txt >> $ofile
        echo "<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\"><PRE>" >> $ofile
        grep -v photometry $file >> $ofile
        echo "</PRE><BR>" >> $ofile
        grep photometry $file | awk '{if($4>0) printf("<A HREF=\"../%s/thumb_%d_%d_.html\" ID=\"mainimage\" TARGET=\"_blank\"><IMG SRC=\"../%s/thumb_%d_%d_.jpg\" WIDTH=\"'$sz0'\" CLASS=\"map\" USEMAP=\"#%s_thumb_%d_%d_.map\"></A>\n",$1,$2,$3,$1,$2,$3,$1,$2,$3)}' >> $ofile
        for pg in `grep photometry $file | awk '{printf("%s/thumb_%d_%d_.html\n",$1,$2,$3)}'`; do
            echo "<map name=\"$pg\">" | sed -e 's/\//_/' -e 's/html/map/g' >> $ofile
            grep "Source \\$nid," $savedir/$pg  | grep area >> $ofile
            echo "</map>" >> $ofile
        done
        echo "<BR>" >> $ofile  
        grep photometry $file | awk '{if($4>0 && $5>0) printf("<A HREF=\"../%s/lc_%d.txt\" TARGET=\"_blank\"><IMG WIDTH=300 SRC=\"../%s/lc_%d.jpg\"></A>\n",$1,'$nid',$1,'$nid')}' >> $ofile
        echo "</BODY></HTML>" >> $ofile
    fi
done
rm candidate_ids.txt new_sources 2>/dev/null

cp redux_colibri_AL/stack_AL.jpg redux_colibri_AL/stack_AL_circles.jpg redux_colibri_AL/ps1_stack.jpg $pdir0

pdir=${savedir}/$pdir0
[ -d $pdir ] && rm -r $pdir
mv $index $pdir0 $savedir
