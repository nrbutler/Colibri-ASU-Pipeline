#!/bin/bash

snr_cut=5.0
thumb_size=401
sz=1024
savedir=`pwd`
purge=no   # purge single catalog detections

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

hdr0="#    id    RA           DEC                             mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname"
hdr="#    id   RA           DEC                                 mag          dmag         fwhm high_pm fade subdet deblend grb_offset catmag catoffset catname"

#
# collect sources from the individual webpages
#

rm bright_fading.txt bright_nonfading.txt bright_questionable.txt 2>/dev/null
nfiles0=`cat redux_colibri_C?_*/fmatched*_list.txt | wc -l`
filters=
for hfile in `ls redux_colibri_C?_*/index.html`; do
    photdir=`grep index_faint $hfile | awk '{print $2}' | sed -e 's/"/ /g' -e 's/\// /g' | awk '{print $2}'`
    filter=`echo $photdir | awk -F_ '{print $(NF-1)}'`
    filters="${filters}$filter"
    awk '/Bright\/Fading/,/Bright\/Non-Fading/{print}' $hfile | egrep -v 'dslope|Fading' | sed -e "s/<A/${filter}<A/" >> bright_fading.txt
    awk '/Bright\/Non-Fading/,/Questionable/{print}' $hfile | egrep -v 'dslope|Bright' | sed -e "s/<A/${filter}<A/" >> bright_nonfading.txt
    awk '/Questionable/,/sourcelists_done/{print}' $hfile | egrep -v 'dslope|Bright|sourcelists' | sed -e "s/<A/${filter}<A/" >> bright_questionable.txt
    dir=`dirname $hfile`
    grep HREF ${dir}/index_faint.html | grep thumb | sed -e "s/thumb/${photdir}\/thumb/g" -e "s/<A/${filter}<A/" -e 's/decoration:none\">/decoration:none\"> /g' >> faint_sources.txt
done

if [ -f faint_sources.txt ]; then
    sort -nr -k 5 faint_sources.txt > faint_sources.tmp
    mv faint_sources.tmp faint_sources.txt
fi
if [ -f bright_fading.txt ]; then
    sort -nr -k 4 bright_fading.txt > bright_fading.tmp
    mv bright_fading.tmp bright_fading.txt
fi
if [ -f bright_nonfading.txt ]; then
    sort -nr -k 4 bright_nonfading.txt > bright_nonfading.tmp
    mv bright_nonfading.tmp bright_nonfading.txt
fi
if [ -f bright_questionable.txt ]; then
    sort -nr -k 4 bright_questionable.txt > bright_questionable.tmp
    mv bright_questionable.tmp bright_questionable.txt
fi

tag="${filters}_${nfiles0}"
pdir0=summary_$tag
rm -r $pdir0 2>/dev/null
rm new_sources 2>/dev/null
mkdir $pdir0
ln -s $pdir0 new_sources

# webpage to make
index=summary_${tag}.html
index_faint=new_sources/faint_${tag}.html
index_cat=new_sources/cat_${tag}.html
index_quest=new_sources/quest_${tag}.html
rm $index $index_faint $index_cat 2>/dev/null

img_circ=`grep 'IMG SRC' redux_colibri_AL/index.html | awk -F\" '{print $2}'`
img_nocirc=`echo $img_circ | sed -e 's/_circles//g'`
ps1file_jpeg=`grep 'PS1 image' redux_colibri_AL/index.html | awk -F\" '{print $2}'`
difffile_jpeg=`grep 'Difference image' redux_colibri_AL/index.html | awk -F\" '{print $2}'`

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
<BUTTON ONCLICK=diffImage(\"${img_nocirc}\")>Image without circles</BUTTON> &nbsp; &nbsp;" > $index
[ "$ps1file_jpeg" ] && echo "<BUTTON ONCLICK=diffImage(\"${ps1file_jpeg}\")>PS1 image</BUTTON> &nbsp; &nbsp;" >> $index
[ "$difffile_jpeg" ] && echo "<BUTTON ONCLICK=diffImage(\"${difffile_jpeg}\")>Difference image</BUTTON> &nbsp; &nbsp;" >> $index

echo "<P>" >> $index
awk '/COLIBRI C/,/FWHM:/{print}' redux_colibri_C?_*/index.html >> $index
awk '/COLIBRI AL/,/FWHM:/{print}' redux_colibri_AL/index.html >> $index

echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">Bright/Fading (SNR>=$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> $index
echo "$hdr0" >> $index
awk '{id=$4; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_fading.txt | sed -e 's/,\.html/\.html/g' >> $index
echo "</PRE><FONT SIZE=\"+1\" COLOR=\"#006600\">Bright/Non-Fading (SNR>=$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" >> $index
echo "$hdr0" >> $index
awk '{id=$4; printf("%s <A HREF=\"'$pdir0'/source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_nonfading.txt | sed -e 's/,\.html/\.html/g' >> $index
echo "</PRE><A HREF=\"${pdir0}/faint_${tag}.html\" TARGET=\"_blank\">Faint New Source Photometry</A>
&nbsp; &nbsp; <A HREF=\"${pdir0}/quest_${tag}.html\" TARGET=\"_blank\">Questionable Source Photometry</A>
&nbsp; &nbsp; <A HREF=\"${pdir0}/det_effic.jpg\" TARGET=\"_blank\"> Detection Efficiency </A><BR>" >> $index

# make a detection efficiency plot
get_det_effic.py >> $index
if [ -f det_effic.jpg ]; then
    mv det_effic.jpg $pdir0
    echo "<BR>" >> $index
fi

if [ "$purge" = "yes" ]; then
    echo "<A HREF=\"${pdir0}/cat_${tag}.html\" TARGET=\"_blank\">Partially Catalogued Source Photometry</A><BR>" >> $index
fi

echo "<HTML><HEAD><TITLE>Questionable Sources</TITLE></HEAD><BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\"
<FONT SIZE=\"+1\" COLOR=\"#006600\">Questionable (SNR>$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" > $index_quest
echo "$hdr0" >> $index_quest
awk '{id=$4; printf("%s <A HREF=\"source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' bright_questionable.txt | sed -e 's/,\.html/\.html/g' -e 's/phot/..\/phot/g' >> $index_quest
echo "</PRE></BODY></HTML>" >> $index_quest

echo "<HTML><HEAD><TITLE>Faint Sources</TITLE></HEAD><BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\"
<FONT SIZE=\"+1\" COLOR=\"#006600\">Faint (SNR<$snr_cut) Uncatalogued Sources:</FONT><BR><P><PRE>" > $index_faint
echo "$hdr" >> $index_faint
awk '{id=$5; printf("%s <A HREF=\"source_%s.html\" TARGET=\"_blank\">[full source info]</A>\n",$0,id)}' faint_sources.txt | sed -e 's/<\/A>\.html/\.html/g' -e 's/phot/..\/phot/g' >> $index_faint
echo "</PRE></BODY></HTML>" >> $index_faint

if [ "$purge" = "yes" ]; then
    echo "<HTML><HEAD><TITLE>Catalogued Sources</TITLE></HEAD><BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\"
        <FONT SIZE=\"+1\" COLOR=\"#006600\">Partially Catalogued Sources:</FONT><BR><P><PRE>" > $index_cat
fi

# retrieve the candidate id numbers
cat bright_fading.txt bright_nonfading.txt | awk '{print $4,0}' | sed -e 's/,//g' > candidate_ids.txt
[ -f bright_questionable.txt ] && cat bright_questionable.txt | awk '{print $4,2}' | sed -e 's/,//g' >> candidate_ids.txt
awk '{print $5,1}' faint_sources.txt | sed -e 's/<\/A>//' >> candidate_ids.txt
# get unique ids, questionable overrides faint, bright overrides all
awk '{id[$1]=$1; if($2==0) q[$1]=0; if(!($1 in q) || q[$1]==1) q[$1]=$2}END{for (indx in id) {print id[indx],q[indx]}}' candidate_ids.txt | sort -nr -k 1 > ids$$.tmp
mv ids$$.tmp candidate_ids.txt

# compile source info
source_info.py candidate_ids.txt $snr_cut

# now purge out sources catalogued in any epoch
if [ "$purge" = "yes" ]; then
    grep Catalog new_sources/source_-*.html | grep -v NONE | awk -F/ '{print $2}' | awk -F: '{print $1}' > purge$$.txt
    awk '{printf("source_%d.txt %d\n",$1,$2)}' candidate_ids.txt > candidate_id_files.txt
    rm index_cat.tmp 2>/dev/null
    for file in `cat purge$$.txt`; do

        grep $file $index >> index_cat.tmp
        grep -v $file $index > index$$.tmp
        mv index$$.tmp $index

        grep $file $index_faint >> index_cat.tmp
        grep -v $file $index_faint > index_faint$$.tmp
        mv index_faint$$.tmp $index_faint

        grep -v $file candidate_id_files.txt > candidate_id_files.tmp
        mv candidate_id_files.tmp candidate_id_files.txt 

    done
    sed -e 's/source_//g' -e 's/\.txt//g' candidate_id_files.txt > candidate_ids.txt
    awk '{print $0 > NF"_cat.tmp"}' index_cat.tmp
    echo "$hdr0" >> $index_cat
    cat 34_cat.tmp 2>/dev/null | sed -e 's/phot/..\/phot/g' -e 's/summary_.*\/source/source/g' >> $index_cat
    echo "</PRE><P><PRE>$hdr" >> $index_cat
    cat 25_cat.tmp 2>/dev/null >> $index_cat
    echo "</PRE></BODY></HTML>" >> $index_cat
    rm purge$$.txt candidate_id_files.txt index_cat.tmp *_cat.tmp 2>/dev/null
fi

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
    rm radec0_$$.txt
fi

xy_bright=xy0_$$.txt
xy_faint=xy1_$$.txt
xy_quest=xy2_$$.txt
#egrep 'Source|Catalog' $pdir0/source_*.html | awk '{if($0~/NONE/) print dat; dat=$0}' > sources$$.txt
grep 'Source' $pdir0/source_*.html | grep -v TITLE > sources$$.txt
grep Bright sources$$.txt | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_bright
grep Faint sources$$.txt | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_faint
grep Questionable sources$$.txt | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_quest
rm sources$$.txt 2>/dev/null

#grep Source $pdir0/source_*.html 2>/dev/null | grep Bright | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_bright
#grep Source $pdir0/source_*.html 2>/dev/null | grep Faint | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_faint
#grep Source $pdir0/source_*.html 2>/dev/null | grep Questionable | awk '{print $7,$2}' | tr '(,:)' ' ' > $xy_quest

redcircles0=
[ -s "$xy_bright" ] && redcircles0=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';n=$3; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -stroke black -fill green -draw '\''translate %.0f,%.0f text -5,-11 \"%.0f\"'\''",x,y,x,y+5,x,y,n)}' $xy_bright`
redcircles=
[ -s "$xy_faint" ] && redcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\''",x,y,x,y+5)}' $xy_faint`
cyancircles=
[ -s "$xy_quest" ] && cyancircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke cyan -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\''",x,y,x,y+5)}' $xy_quest`

echo "<map name=\"summary.map\">" > summary.map
cat $xy_bright $xy_faint $xy_quest | awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$3; printf("<area shape =\"circle\" coords=\"%.0f,%.0f,5\" href=\"'${pdir0}'/source_%.0f.html\" target=\"_blank\" title=\"Source %.0f\">\n",x,y,id,id)}' >> summary.map
echo "</map>" >> summary.map
rm $xy_bright $xy_faint $xy_quest 2>/dev/null

if [ -f summary.map ]; then
    cat summary.map >> $index
    rm summary.map
fi

# make a detection efficiency plot
get_det_effic.py
[ -f det_effic.jpg ] && mv det_effic.jpg $pdir0

echo "<P>Previous Reductions: " >> $index
for redux in `ls $savedir/summary_*html 2>/dev/null | grep -v $index`; do
    redux=`basename $redux`
    num=`echo $redux | awk -F_ '{print $3}' | sed -e 's/\.html//g'`
    echo "[<A HREF=\"$redux\" TARGET=\"_blank\"> $num </A>]" 
done | sort -nr -k 4 >> $index
dte=`date -u`

here=`pwd`
cd $savedir
rm listing.html 2>/dev/null
make_listing.sh | grep -v listing > listing.html
cd $here

echo "<P><A HREF=\"listing.html\" TARGET=\"_blank\">Directory Listing</A>
<P><HR WIDTH=\"100%\">Last Updated: $dte (natbutler@asu.edu)</BODY></HTML>" >> $index

cmd="convert -quality 75 redux_colibri_AL/stack_AL.jpg -pointsize 16 $redcircles0 $redcircles $cyancircles $greencircles ${pdir0}/stack_AL_circles.jpg"
eval "$cmd"

# boldface repeat entries in html files
for html_file in `ls $index $index_faint $index_cat $index_quest 2>/dev/null`; do
    grep 'full source info' $html_file | awk '{print $(NF-3)}' | sort | uniq -c | awk '{if($1>1) print $2}' | sed -e 's/"//g' -e 's/\//\\\//g' -e 's/HREF=//g' > repeat_tags.txt
    for tag in `cat repeat_tags.txt`; do
        awk '{if($0~/'$tag'/ && $0~/A HREF/) print "<B>"$0"</B>"; else print}' $html_file > html_file.tmp
        mv html_file.tmp $html_file
    done
    rm repeat_tags.txt 2>/dev/null
done

pdir=${savedir}/$pdir0
[ -d $pdir ] && rm -r $pdir
mv $index $pdir0 $savedir

rm faint_sources.txt bright_fading.txt bright_nonfading.txt bright_questionable.txt candidate_ids.txt new_sources 2>/dev/null

cd $savedir
rm index.html 2>/dev/null
ln -s $index index.html
