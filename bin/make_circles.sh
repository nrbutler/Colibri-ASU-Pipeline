#!/bin/sh

fitsfile=$1
shift
radeclist=$1
shift

sz=1024
scl=1

x0=`gethead CRPIX1 $fitsfile`
y0=`gethead CRPIX2 $fitsfile`
calfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

x=`gethead CRPIX1 $fitsfile`
y=`gethead CRPIX2 $fitsfile`

nx=`gethead NAXIS1 $fitsfile`
ny=`gethead NAXIS2 $fitsfile`
sx=$sz
sy=`echo "$sz $nx $ny" | awk '{printf("%.0f\n",$1*$3/$2)}'`

facx=`echo $sx $nx | awk '{printf("%.6f\n",$2/$1)}'`
facy=`echo $sy $ny | awk '{printf("%.6f\n",$2/$1)}'`

pos_err=`gethead SEXERR $fitsfile`
if [ "$pos_err" ]; then
    xpos=`gethead SEXX0 $fitsfile`
    ypos=`gethead SEXY0 $fitsfile`
    echo $xpos $ypos $pos_err > radec0_$$.txt
fi

base=${fitsfile%'.fits'}
[ -f ${base}.jpg ] || fits2jpg.py $fitsfile $sz linvert

echo "<map name=\"${base}_circles.map\">" > ${base}_circles.map
echo "<map name=\"${base}.map\">" > ${base}.map


if [ -f "$radeclist" ]; then
    radeclist1=f$radeclist
    grep -v '#' $radeclist > $radeclist1
    if [ -s "$calfile" ]; then
        calfile1=f$calfile
        grep -v '#' $calfile | awk '{if($NF<0) print}' >> $radeclist1
        grep -v '#' $calfile | awk '{if($NF>=0) print}' > $calfile1
    fi
fi

redcircles=
if [ -f "$radeclist1" ]; then
    awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,$3,$NF}' $radeclist1 > radec$$.txt

    redcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';n=$NF; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\''",x,y,x,y+1)}' radec$$.txt`

    nf=`cat radec$$.txt | wc -l`
    ls lc_dir/lc_*txt | sed -e 's/lc_dir\/lc_//g' -e 's/\.txt//g' | sort -nr |\
       awk 'BEGIN{n=-1}{for(i=n;i>$1;i--) print ")\">"; printf(", click for lc)\" href=\"lc_%d.jpg\">\n",$1); n=$1-1}END{for(i=n;i>=-'$nf';i--) print ")\">"}' > lc_link$$.txt
    awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$NF; printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" target=\"_blank\" title=\"Source %d, mag=%.1f (%.0f %.0f\n",x*'$scl',y*'$scl',id,$3,$1,$2)}' radec$$.txt > ${base}.tmp
    paste -d "" ${base}.tmp lc_link$$.txt >> ${base}.tmp1
    cat ${base}.tmp1 >> ${base}_circles.map
    cat ${base}.tmp1 >> ${base}.map

    #awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$NF; if (system("[ -f lc_dir/lc_"id".jpg ]")) printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" target=\"_blank\" title=\"Source %d, mag=%.1f (%.0f %.0f)\">\n",x*'$scl',y*'$scl',id,$3,$1,$2); else printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" href=\"lc_%d.jpg\" target=\"_blank\" title=\"Source %d, mag=%.1f (%.0f %.0f, click for lc)\">\n",x*'$scl',y*'$scl',id,id,$3,$1,$2)}' radec$$.txt > ${base}.tmp
    #cat ${base}.tmp >> ${base}_circles.map
    #cat ${base}.tmp >> ${base}.map

    rm radec$$.txt ${base}.tmp ${base}.tmp1 2>/dev/null
fi

bluecircles=
if [ "$calfile1" ]; then
    if [ -s "$calfile1" ]; then
        awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,$3,$NF}' $calfile1 > radec$$.txt
        bluecircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke blue -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' ",x,y,x,y+2)}' radec$$.txt`
        #awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$NF; if (system("[ -f lc_dir/lc_"id".jpg ]")) printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" target=\"_blank\" title=\"USNO %d, mag=%.1f (%.0f %.0f)\">\n",x*'$scl',y*'$scl',id,$3,$1,$2); else printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" href=\"lc_%d.jpg\" target=\"_blank\" title=\"USNO %d, mag=%.1f (%.0f %.0f, click for lc)\">\n",x*'$scl',y*'$scl',id,id,$3,$1,$2)}' radec$$.txt > ${base}.tmp
        awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; id=$NF; printf("<area shape=\"circle\" coords=\"%.0f,%.0f,5\" target=\"_blank\" title=\"USNO %d, mag=%.1f (%.0f %.0f)\">\n",x*'$scl',y*'$scl',id,$3,$1,$2)}' radec$$.txt > ${base}.tmp
        cat ${base}.tmp >> ${base}_circles.map
        cat ${base}.tmp >> ${base}.map
        rm radec$$.txt ${base}.tmp 2>/dev/null
    fi
fi

yellowcircles=
if [ "$pos_err" ]; then
    yellowcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke yellow -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -draw '\''line 0,%.0f %.0f,%.0f'\'' -draw '\''line %.0f,0 %.0f,%.0f'\'' ",x,y,x,y+$3/'$facy',y,('$nx'-1)/'$facx',y,x,x,('$ny'-1)/'$facy')}' radec0_$$.txt`
    rm radec0_$$.txt
fi

sz0=`echo $sz $scl | awk '{print $1*$2}'`
echo "<area shape=\"rect\" coords=\"0,0,$sz0,$sz0\" href=\"${base}_circles.html\" title=\"Click to Return Circles\"></map>" >> ${base}.map
echo "<area shape=\"rect\" coords=\"0,0,$sz0,$sz0\" href=\"${base}.html\" title=\"Click to Remove Circles\"></map>" >> ${base}_circles.map

cmd="convert -quality 75 ${base}.jpg -pointsize 16 $redcircles $bluecircles $yellowcircles ${base}_circles.jpg"
eval "$cmd"

# remove red circles (new sources) if it was unable to make them
cmd="convert -quality 75 ${base}.jpg -pointsize 16 $bluecircles $yellowcircles ${base}_circles.jpg"
[ -f "${base}_circles.jpg" ] || eval "$cmd"

rm $radeclist1 $calfile1
