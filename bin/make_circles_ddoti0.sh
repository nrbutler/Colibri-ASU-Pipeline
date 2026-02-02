#!/bin/bash

fitsfile=$1
shift
radeclist=$1
shift
radeclist_faint=$1
shift
radeclist_catfaint=$1
shift
radeclist_all=$1
shift

sz=1024
thumb_size=401
thumb_edge=30
thumb_scale=1
i0=
j0=

make_thumbs=yes
make_fits_thumbs=no

logscale=yes
invert_cmap=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

nx=`gethead NAXIS1 $fitsfile`
ny=`gethead NAXIS2 $fitsfile`
sx=$sz
sy=`echo "$sz $nx $ny" | awk '{printf("%.0f\n",$1*$3/$2)}'`

filter=`gethead FILTER $fitsfile`
[ "$filter" ] || filter=r
base=${fitsfile%'.fits'}

diff_file=${base}.diff.fits

pfilter=$filter
[ "$pfilter" = "B" ] && pfilter=g
[ "$pfilter" = "gri" ] && pfilter=r
[ "$pfilter" = "zy" ] && pfilter=z
ps1_file=ps1_stack_${pfilter}.fits

facx=`echo $sx $nx | awk '{printf("%.6f\n",$2/$1)}'`
facy=`echo $sy $ny | awk '{printf("%.6f\n",$2/$1)}'`

# defined in fits2jpg.py:
[ "$i0" ] || i0=`echo $nx $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`
[ "$j0" ] || j0=`echo $ny $thumb_size | awk '{n=int($1/$2); if(n<$1/$2) n=n+1; print int(($1-n*$2)/2)}'`

fov=`echo $thumb_size $thumb_edge | awk '{printf("%.6f\n",($1+2*$2)*0.38/3600)}'`

#TNS
turl1="https://www.wis-tns.org/search?ra="
turl2="&radius=10&coords_unit=arcsec&include_frb=1"

xpos=-1; ypos=-1; xthumb=-1; ythumb=-1
pos_err=`gethead SEXERR $fitsfile`
if [ "$pos_err" ]; then
    xpos=`gethead SEXX0 $fitsfile`
    ypos=`gethead SEXY0 $fitsfile`
    xthumb=`echo "$xpos" | awk 'BEGIN{t='$thumb_size'}{print '$i0'+int(($1-('$i0'))/t)*t}'`
    ythumb=`echo "$ypos" | awk 'BEGIN{t='$thumb_size'}{print '$j0'+int(($1-('$j0'))/t)*t}'`
    echo $xpos $ypos $pos_err > radec0_$$.txt
fi

inv=noinvert
[ "$invert_cmap" = "yes" ] && inv=invert
[ "$logscale" = "yes" ] && inv=l$inv

mt=
if [ "$make_thumbs" = "yes" ]; then
    rm thumb_*_*_*_*_*_.jpg 2>/dev/null
    mt=make_thumbs
fi

# create jpgs (and possibly fits thumbnails)
save_fits=no
[ "$make_fits_thumbs" = "yes" ] && save_fits=save_fits
echo fits2jpg_ddoti.py $fitsfile $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits $filter
fits2jpg_ddoti.py $fitsfile $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits $filter

if [ -f "$diff_file" ]; then
    echo fits2jpg_ddoti.py $diff_file $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits $filter diff_
    fits2jpg_ddoti.py $diff_file $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits $filter diff_
else
    echo "No image subtraction file $diff_file found."
fi

if [ -f "$ps1_file" ]; then
    echo fits2jpg_ddoti.py $ps1_file $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits PS1$pfilter ps1_
    fits2jpg_ddoti.py $ps1_file $sz $inv $mt $thumb_size $thumb_edge [$j0,$i0] $save_fits PS1$pfilter ps1_
else
    echo "No PS1 file $ps1_file found."
fi

sz0=`echo $thumb_size $thumb_edge $thumb_scale | awk '{print ($1+2*$2)*$3}'`
C11=`gethead CD1_1 $fitsfile`; C12=`gethead CD1_2 $fitsfile`
pscale=`echo $C11 $C12 | awk '{print sqrt($1*$1+$2*$2)*3600.}'`
phys_size=`echo $thumb_size $thumb_edge $pscale | awk '{printf("%.1f\n",($1+2*$2)*$3/60.)}'`
phys_edge=`echo $thumb_edge $pscale | awk '{printf("%.1f\n",$1*$2/60.)}'`

#
# mark up the thumbs
#
if [ "$make_thumbs" = "yes" ]; then

    # make an html page for each thumb jpg
    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{print "<HTML>">"thumb_"$2"_"$3"_.html"}'

    for file in `ls thumb_*_*_.html`; do
        cat ${REDUX_BASE_DIR}/jquery/hdr.txt >> $file
    done

    # make some symbolic links
    ls thumb_*_*_*_*_*_.jpg | awk -F_ '{printf("ln -s %s thumb_%d_%d_.jpg\n",$0,$2,$3)}' | sh
    ls ps1_thumb_*_*_*_*_*_.jpg 2>/dev/null | awk -F_ '{printf("ln -s %s ps1_thumb_%d_%d_.jpg\n",$0,$3,$4)}' | sh
    ls diff_thumb_*_*_*_*_*_.jpg 2>/dev/null | awk -F_ '{printf("ln -s %s diff_thumb_%d_%d_.jpg\n",$0,$3,$4)}' | sh

    # build image maps
    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("<BODY>\n<map name=\"thumb_%d_%d_.map\">\n",$2,$3)>>"thumb_"$2"_"$3"_.html"}'
    grep -v '#' $radeclist_all 2>/dev/null | sort -n -k 5 | awk 'BEGIN{t='$thumb_size';ts='$thumb_scale';te='$thumb_edge'}{color="0000ff"; if ($7==1) {color="aa0000"; if ($8==1) color="ff0000"}; rad=$6*ts; if(rad<2*ts) rad=2*ts ;x='$i0'+int(($1-('$i0')-0.5)/t)*t; y='$j0'+int(($2-('$j0'))/t)*t; if ($8==0) printf("<area shape=\"circle\" coords=\"%.0f,%.0f,%.0f\" href=\"#\" onclick=\"return false;\" title=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" id=\"%d\" data-maphilight=\047{\"strokeColor\":\"%s\",\"strokeWidth\":\"2\",\"alwaysOn\":\"true\",\"fillOpacity\":\"0\"}\047>\n",($1-x-1+te)*ts,(t+te-$2+y)*ts,rad,$5,$3,$4,$10,$11,$9,$5,color)>>"thumb_"x"_"y"_.html"; else printf("<area shape=\"circle\" coords=\"%.0f,%.0f,%.0f\" href=\"lc_%d.jpg\" target=\"_blank\" title=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" id=\"%d\" data-maphilight=\047{\"strokeColor\":\"%s\",\"strokeWidth\":\"2\",\"alwaysOn\":\"true\",\"fillOpacity\":\"0\"}\047>\n",($1-x-1+te)*ts,(t+te-$2+y)*ts,rad,$5,$5,$3,$4,$10,$11,$9,$5,color)>>"thumb_"x"_"y"_.html"}'

    # make the source lists
    rm radec_*_*_.txt 2>/dev/null
    grep -v '#' $radeclist 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{x='$i0'+int(($6-('$i0')-0.5)/t)*t; y='$j0'+int(($7-('$j0')-0.5)/t)*t; print $6-x,$7-y,$8,$1>"radec_"x"_"y"_.txt"}'

    function circle_thumb() {
        local base=$1
        local radecfile=$2
        local id
        local x0=`echo $base | awk -F_ '{print $2}'`
        local y0=`echo $base | awk -F_ '{print $3}'`
        local base0=`echo $base | awk -F_ '{print $1"_"$2"_"$3"_"}'`

        local redcircles=
        if [ -f "$radecfile" ]; then
          local redcircles=`awk '{x=$1-1+'$thumb_edge';y='$thumb_size'+'$thumb_edge'-$2;n=$4;\
            printf(" -stroke black -fill green -draw '\''translate %.0f,%.0f text -5,-11 \"%d\"'\''",x,y,n)}' $radecfile`
        fi
 
        local greencircles=
        xin=0; yin=0
        if [ "$pos_err" ]; then
            xin=`echo $thumb_size $thumb_edge $pos_err $xpos $x0 | awk '{xin=0;ts=$1;te=$2;pe=$3;dx=$4-$5;if(dx+pe>=1-te && dx-pe<=ts+te) xin=1; print xin}'`
            yin=`echo $thumb_size $thumb_edge $pos_err $ypos $y0 | awk '{yin=0;ts=$1;te=$2;pe=$3;dy=$5-$4;if(dy+pe>=-ts-te && dy-pe<=te-1) yin=1; print yin}'`
        fi
        #if [ "$x0" -eq "$xthumb" -a "$y0" -eq "$ythumb" ]; then
        if [ "$xin" -eq 1 -o "$yin" -eq 1 ]; then
            greencircles=`awk '{x=($1-'$x0'-1+'$thumb_edge');y=('$thumb_size'+'$thumb_edge'-$2+'$y0'); printf(" -stroke green -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -draw '\''line 0,%.0f %.0f,%.0f'\'' -draw '\''line %.0f,0 %.0f,%.0f'\'' ",x,y,x,y+$3,y,('$thumb_size'+2*'$thumb_edge'-1),y,x,x,('$thumb_size'+2*'$thumb_edge'-1))}' radec0_$$.txt`
        fi

        local cmd="convert -quality 75 ${base}.jpg -pointsize 16 $redcircles $greencircles ${base}_circles.jpg"
        eval "$cmd"
        mv ${base}.jpg ${base}nocircles.jpg
        if [ -f "${base}_circles.jpg" ]; then
            mv ${base}_circles.jpg ${base}.jpg
        else
            cp ${base}nocircles.jpg ${base}.jpg
        fi
    }

    if [ "${xthumb}" -ge 0 -a "${ythumb}" -ge 0 ]; then
        files="radec_${xthumb}_${ythumb}_.txt "`ls radec_*_*_.txt 2>/dev/null | grep -v "radec_${xthumb}_${ythumb}_.txt"`
    else
        files=`ls radec_*_*_.txt 2>/dev/null`
    fi

    for radecfile in $files; do
        x=`echo $radecfile | awk -F_ '{print $2}'`
        y=`echo $radecfile | awk -F_ '{print $3}'`
        jfile=`ls thumb_${x}_${y}_*_*_*_.jpg`
        r=`echo $jfile | awk -F_ '{print $4}'`
        d=`echo $jfile | awk -F_ '{print $5}'`
        th=`echo $jfile | awk -F_ '{print $6}'`
        circle_thumb thumb_${x}_${y}_${r}_${d}_${th}_ radec_${x}_${y}_.txt &
    done
    wait

    # finish off each webpage
    if [ -f "$ps1_file" ]; then
        ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{r=$4;d=$5;gsub(":"," ",r);gsub(":"," ",d); radec=r" "d; printf("</map>\n <table BORDER=\"0\" CELLPADDING=\"32\" CELLSPACING=\"0\" BGCOLOR=\"#ffffff\">\n<tr><td>\n<IMG SRC=\"%s.jpg\" ID=\"mainimage\" USEMAP=\"#thumb_%d_%d_.map\" WIDTH=\"'$sz0'\" CLASS=\"map\"> \n </td><td>\n <IMG SRC=\"ps1_%s.jpg\" ID=\"mainimage\" USEMAP=\"#thumb_%d_%d_.map\" WIDTH=\"'$sz0'\" CLASS=\"map\"> \n </td></tr>\n </table>\n <BR> <BUTTON ONCLICK=change_table()>Switch Images</Button> &nbsp;&nbsp;&nbsp; Physical Size: '$phys_size' arcmin. &nbsp;&nbsp;&nbsp; <A ID=\"toggleall\" HREF=\"#\">Toggle All Sources</A> &nbsp;&nbsp;&nbsp; <A HREF=\"thumb_"($2-'$thumb_size')"_"$3"_.html\" target=\"_blank\">Thumb Left</A> <A HREF=\"thumb_"($2+'$thumb_size')"_"$3"_.html\" target=\"_blank\">Thumb Right</A> <A HREF=\"thumb_"$2"_"($3+'$thumb_size')"_.html\" target=\"_blank\">Thumb Up</A> <A HREF=\"thumb_"$2"_"($3-'$thumb_size')"_.html\" target=\"_blank\">Thumb Down</A> <P> Note: sources in outer '$phys_edge' arcmin of image are displayed on the neighboring thumbnail. &nbsp;&nbsp;&nbsp; <A HREF=\"thumb_"$2"_"$3"_.html\" target=\"_blank\">Unmarked Thumb</A> &nbsp;&nbsp;&nbsp; <A HREF=\"thumb_"$2"_"$3"_.html\" target=\"_blank\">Difference Image Thumb</A> <BR><HR>Bright Uncatalogued Sources:<BR><PRE>\n#   id        RA           DEC                              mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname\n",$0,$2,$3,$0,$2,$3)>>"thumb_"$2"_"$3"_.html"}'
    else
        ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{r=$4;d=$5;gsub(":"," ",r);gsub(":"," ",d); radec=r" "d; printf("</map>\n <table BORDER=\"0\" CELLPADDING=\"32\" CELLSPACING=\"0\" BGCOLOR=\"#ffffff\">\n<tr><td>\n<IMG SRC=\"%s.jpg\" ID=\"mainimage\" USEMAP=\"#thumb_%d_%d_.map\" WIDTH=\"'$sz0'\" CLASS=\"map\"> \n </td><td>\n <div id=\"aladin-lite-div\" style=\"width:'$sz0'px;height:'$sz0'px;\"></div> <script type=\"text/javascript\" src=\"http://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js\" charset=\"utf-8\"></script> <script type=\"text/javascript\"> var aladin = A.aladin(\"#aladin-lite-div\", {target: \"%s\", survey: \"P/PanSTARRS/DR1/g\", fov: '$fov', showReticle: false}); aladin.getBaseImageLayer().getColorMap().reverse(); </script>\n <script type=\"text/javascript\"> $(\"#aladin-lite-div\").rotate(%d); </script>\n </td></tr>\n </table>\n <BR> <BUTTON ONCLICK=change_table()>Switch Images</Button> &nbsp;&nbsp;&nbsp; Physical Size: '$phys_size' arcmin. &nbsp;&nbsp;&nbsp; <A ID=\"toggleall\" HREF=\"#\">Toggle All Sources</A> &nbsp;&nbsp;&nbsp; <A HREF=\"thumb_"($2-'$thumb_size')"_"$3"_.html\" target=\"_blank\">Thumb Left</A> <A HREF=\"thumb_"($2+'$thumb_size')"_"$3"_.html\" target=\"_blank\">Thumb Right</A> <A HREF=\"thumb_"$2"_"($3+'$thumb_size')"_.html\" target=\"_blank\">Thumb Up</A> <A HREF=\"thumb_"$2"_"($3-'$thumb_size')"_.html\" target=\"_blank\">Thumb Down</A> <P> Note: sources in outer '$phys_edge' arcmin of image are displayed on the neighboring thumbnail. &nbsp;&nbsp;&nbsp; <A HREF=\"thumb_"$2"_"$3"_.html\" target=\"_blank\">Unmarked Thumb</A> <BR><HR>Bright Uncatalogued Sources:<BR><PRE>\n#   id        RA           DEC                              mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname\n",$0,$2,$3,radec,$6)>>"thumb_"$2"_"$3"_.html"}'
    fi

    fwhm0=`gethead FWHM $fitsfile`

    grep -v '#' $radeclist 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{good=1; if($10>3*'$fwhm0') good=0; if ($13>10 && sqrt($11*$11)<0.1) good=0; if ($8<13) good=0; split($14,nepoch,"/"); if (nepoch[1]<0.75*nepoch[2] && $11>-$12) good=0; if ($15!=0 || $20>0) good=0; if(good==1) {x='$i0'+int(($6-('$i0')-0.5)/t)*t;y='$j0'+int(($7-('$j0')-0.5)/t)*t; r=$2;d=$3;gsub(":","%3A",r); gsub(":","%3A",d); gsub("\+","%2B",d); url="'$turl1'"r"&decl="d"'$turl2'"; printf("<A ID=\"%d\" HREF=\"#\" TITLE=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" TARGET=\"_blank\" style=\"text-decoration:none\">%6d</A> <A HREF=\"lc_%d.jpg\" TARGET=\"_blank\" style=\"text-decoration:none\"> lc </A> %11s %11s (%9.5f,%9.5f) %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %9s %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s <A HREF=\"%s\" target=\"_blank\" TITLE=\"Search TNS\">TNS</A>\n",$1,$1,$8,$9,$21,$22,$20,$1,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,url)>>"thumb_"x"_"y"_.html"}}'

    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("</PRE>Questionable Sources:<BR><PRE>\n#   id        RA           DEC                              mag      dmag     fwhm     slope    dslope       chi2 n_detect high_pm fade subdet deblend grb_offset catmag catoffset catname\n")>>"thumb_"$2"_"$3"_.html"}'

    grep -v '#' $radeclist 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{good=1; if($10>3*'$fwhm0') good=0; if ($13>10 && sqrt($11*$11)<0.1) good=0; if ($8<13) good=0; split($14,nepoch,"/"); if (nepoch[1]<0.75*nepoch[2] && $11>-$12) good=0; if ($15!=0 || $20>0) good=0; if(good==0) {x='$i0'+int(($6-('$i0')-0.5)/t)*t;y='$j0'+int(($7-('$j0')-0.5)/t)*t; r=$2;d=$3;gsub(":","%3A",r); gsub(":","%3A",d); gsub("\+","%2B",d); url="'$turl1'"r"&decl="d"'$turl2'"; printf("<A ID=\"%d\" HREF=\"#\" TITLE=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" TARGET=\"_blank\" style=\"text-decoration:none\">%6d</A> <A HREF=\"lc_%d.jpg\" TARGET=\"_blank\" style=\"text-decoration:none\"> lc </A> %11s %11s (%9.5f,%9.5f) %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %9s %5d %5d %5.2f %5d %8.1f %8.1f %8.1f %10s <A HREF=\"%s\" target=\"_blank\" TITLE=\"Search TNS\">TNS</A>\n",$1,$1,$8,$9,$21,$22,$20,$1,$1,$2,$3,$4,$5,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,url)>>"thumb_"x"_"y"_.html"}}'

    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("</PRE>Faint Uncatalogued Sources:<BR><PRE>\n#   id   RA           DEC                              mag      dmag     fwhm  high_pm subdet deblend grb_offset catmag catoffset catname\n")>>"thumb_"$2"_"$3"_.html"}'

    grep -v '#' $radeclist_faint 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{x='$i0'+int(($6-('$i0')-0.5)/t)*t;y='$j0'+int(($7-('$j0')-0.5)/t)*t; r=$2;d=$3; printf("<A ID=\"%d\" HREF=\"#\" TITLE=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" TARGET=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %8.4f %8.4f %8.4f %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",$1,$1,$8,$9,$17,$18,$16,$1,$2,$3,$4,$5,$8,$9,$10,$11,$13,$14,$15,$16,$17,$18)>>"thumb_"x"_"y"_.html"}'

    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("</PRE>Possibly Catalogued Sources:<BR><PRE>\n#   id   RA           DEC                              mag      dmag     fwhm  high_pm subdet deblend grb_offset catmag catoffset catname\n")>>"thumb_"$2"_"$3"_.html"}'

    grep -v '#' $radeclist_catfaint 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{x='$i0'+int(($6-('$i0')-0.5)/t)*t;y='$j0'+int(($7-('$j0')-0.5)/t)*t; r=$2;d=$3; printf("<A ID=\"%d\" HREF=\"#\" TITLE=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" TARGET=\"_blank\" style=\"text-decoration:none\">%6d</A> %11s %11s (%9.5f,%9.5f) %8.4f %8.4f %8.4f %5d %5.2f %5d %8.1f %8.1f %8.1f %10s\n",$1,$1,$8,$9,$17,$18,$16,$1,$2,$3,$4,$5,$8,$9,$10,$11,$13,$14,$15,$16,$17,$18)>>"thumb_"x"_"y"_.html"}'

    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("</PRE><BR>Catalogued Sources:<BR>\n")>>"thumb_"$2"_"$3"_.html"}'

    grep -v '#' $radeclist_all 2>/dev/null | awk 'BEGIN{t='$thumb_size'}{if($7==0) {x='$i0'+int(($1-('$i0')-0.5)/t)*t; y='$j0'+int(($2-('$j0')-0.5)/t)*t; if($5>0) printf("<A ID=\"%d\" HREF=\"#\" TITLE=\"Source %d, m=%.2f&pm;%.2f (+%.1f&rdquo; %s=%.1f)\" STYLE=\"text-decoration:none\">%d</a>\n",$5,$5,$3,$4,$10,$11,$9,$5)>>"thumb_"x"_"y"_.html"}}'

    ls thumb_*_*_*_*_*_.jpg | sed -e 's/\.jpg//' | awk -F_ '{printf("<HR></BODY></HTML>\n")>>"thumb_"$2"_"$3"_.html"}'

    # make a copy for difference image thumbnails
    if [ -f "$ps1_file" ]; then
       for file in `ls thumb_*_*_*_*_*_.jpg 2>/dev/null`; do
           #file1=${file%'.jpg'}nocircles.jpg
           hfile0=`echo "$file" | awk -F_ '{printf("thumb_%s_%s_.html\n",$2,$3)}'`
           hfile1=`echo "$file" | awk -F_ '{printf("diff_thumb_%s_%s_.html\n",$2,$3)}'`
           sed -e '/SRC=\"thumb/s/thumb/diff_thumb/' -e 's/Difference/Regular/' $hfile0 > $hfile1
           sed -e "s/$hfile0/$hfile1/2" $hfile0 > ${hfile0}.tmp
           mv ${hfile0}.tmp $hfile0
       done
    fi

    # make a copy for unmarked thumbnails
    for file in `ls thumb_*_*_*_*_*_nocircles.jpg 2>/dev/null`; do
        hfile0=`echo "$file" | awk -F_ '{printf("thumb_%s_%s_.html\n",$2,$3)}'`
        hfile=`echo "$file" | awk -F_ '{printf("thumb_%s_%s_nocircles.html\n",$2,$3)}'`
        sed -e '/SRC=\"thumb/s/_\.jpg/_nocircles\.jpg/' -e 's/Unmarked/Marked/' $hfile0 > $hfile
        sed -e "s/$hfile0/$hfile/" $hfile0 > ${hfile0}.tmp
        mv ${hfile0}.tmp $hfile0
    done

fi

#
# mark up the main file
#
base0=${fitsfile%'.fits'}

echo "<map name=\"${base0}.map\">" > ${base0}.map
ls thumb_*_*_*_*_*_.jpg | awk -F_ '{t='$thumb_size';fx='$facx';fy='$facy'; printf("<area shape=\"rect\" coords=\"%f,%f,%f,%f\" href=\"thumb_%d_%d_.html\" title=\"%s,%s zoom\" target=\"_blank\">\n",$2/fx,('$ny'-$3-1-t)/fy,($2+t)/fx,('$ny'-$3-1)/fy,$2,$3,$4,$5)}' >> ${base0}.map
echo "</map>" >> ${base0}.map

redcircles=
if [ -f "$radeclist" ]; then
    grep -v '#' $radeclist | awk '{x=$6;y=$7; if (x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,$8,$1}' > radec_${base0}_$$.txt
    redcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';n=$4; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -stroke black -fill green -draw '\''translate %.0f,%.0f text -5,-11 \"%d\"'\''",x,y,x,y+5,x,y,n)}' radec_${base0}_$$.txt`
    rm radec_${base0}_$$.txt
fi

greencircles=
if [ "$pos_err" ]; then
    greencircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke green -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -draw '\''line 0,%.0f %.0f,%.0f'\'' -draw '\''line %.0f,0 %.0f,%.0f'\'' ",x,y,x,y+$3/'$facy',y,('$nx'-1)/'$facx',y,x,x,('$ny'-1)/'$facy')}' radec0_$$.txt`
fi

cmd="convert -quality 75 ${base0}.jpg -pointsize 16 $redcircles $greencircles ${base0}_circles.jpg"
eval "$cmd"

# tag the ps1 thumbnails if present
#for file in `ls ps1_thumb*jpg`; do
#    ( convert -quality 75 $file -stroke green -pointsize 18 -draw "text 10,20 'PS1'" ${file}.tmp ; mv ${file}.tmp $file ; ) &
#done
#wait
