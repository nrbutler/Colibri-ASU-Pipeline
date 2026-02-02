#!/bin/bash

infile=$1
shift

catdir=.
outfile=ps1.fits
setup_only=no
do_wait=no

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f "$infile" ] || { echo "Input file $infile does not exist." ; exit 1 ; }

ps1_files=`get_ps1imagelist.py 0 0 1 f $infile`
echo "ps1_files: $ps1_files"

[ "$ps1_files" ] || { echo "Cannot find ps1 files to download." ; exit 2 ; }

filter=`gethead FILTER $infile`
pdir=${catdir}/ps1_${filter}
[ -d "$pdir" ] || mkdir $pdir

# download the files

url="https://ps1images.stsci.edu/rings.v3.skycell/"
rm ps1_downloading$$.txt 2>/dev/null
wg=
[ "$do_wait" = "yes" ] || wg="-b"
for file in $ps1_files; do
    cell=`echo $file | awk -F. '{print $4}'`
    scell=`echo $file | awk -F. '{print $5}'`
    url1=${url}${cell}/${scell}

    ofile=${pdir}/${file}.fz
    if [ -f "$ofile" ]; then
        echo "$ofile exists from previous run."
    else
        echo "$ofile" >> ps1_downloading$$.txt
        wget "${url1}/$file" $wg -O $ofile >/dev/null 2>&1 &
    fi
    base=${file%'.fits'}
    ofile=${pdir}/${base}.var.fits.fz
    if [ -f "$ofile" ]; then
        echo "$ofile exists from previous run."
    else
        echo "$ofile" >> ps1_downloading$$.txt
        wget "${url1}/${base}.wt.fits" $wg -O $ofile >/dev/null 2>&1 &
    fi
done
if [ "$do_wait" = "yes" ]; then
    wait
else
    for file in `ls wget-log* 2>/dev/null`; do
        done=`grep saved $file`
        if [ "$done" ]; then
            rm $file
        else
            echo $file >> ps1_downloading$$.txt
        fi
    done
    [ -s "ps1_downloading$$.txt" ] && { echo "Not waiting for ps1 stack download, exiting." ; exit 3 ; }
fi
rm ps1_downloading$$.txt 2>/dev/null

rm ps1_list$$.txt 2>/dev/null
for file in $ps1_files; do
    base=${file%'.fits'}
    ifile=${pdir}/${base}.var.fits.fz
    ofile=${pdir}/${base}.var.fits
    [ -f "$ofile" ] || ps1scale.py $ifile &
    ifile=${pdir}/${file}.fz
    ofile=${pdir}/$file
    [ -f "$ofile" ] || ps1scale.py $ifile &
    wait
    wfile=${pdir}/${base}.wt.fits
    [ -f "$wfile" ] || var2weight.py ${pdir}/${base}.var.fits $wfile
    echo "$ofile" >> ps1_list$$.txt
done

if [ "$setup_only" = "yes" ]; then
    echo "Setup performed, exiting."
    rm ps1_list$$.txt 2>/dev/null
    exit 4
fi

# make a stack
if [ -s "ps1_list$$.txt" ]; then

    boutfile=${outfile%'.fits'}
    imhead -f $infile > ${boutfile}.head
    outwfile=${boutfile}.wt.fits

    sargs="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE Y -SUBTRACT_BACK N -WEIGHT_SUFFIX .wt.fits -WEIGHT_TYPE MAP_WEIGHT -BLANK_BADPIXELS N -GAIN_DEFAULT 1.e9"
    echo swarp @ps1_list$$.txt $sargs -COMBINE_TYPE AVERAGE -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME $outwfile
    swarp @ps1_list$$.txt $sargs -COMBINE_TYPE AVERAGE -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME $outwfile
    rm ${boutfile}.head

    [ -f $outfile ] && fits2jpg.py $outfile 1024 linvert

    if [ -f psf_grid.fits ]; then
        npsf=`gethead NBINX psf_grid.fits`
        x0=`gethead CRPIX1 $infile`; y0=`gethead CRPIX2 $infile`
        binfile=${infile%'.fits'}
        [ -d "${boutfile}_dir" ] && rm -r ${boutfile}_dir
        mkdir ${boutfile}_dir
        cp ${binfile}_dir/catalog.fits ${boutfile}_dir
        [ -f psf_grid0.fits ] && rm psf_grid0.fits 
        echo calc_psf $outfile ${boutfile}.wt.fits ${binfile}_dir/seg.fits ${boutfile}_dir/catalog.fits ${binfile}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 0
        calc_psf $outfile ${boutfile}.wt.fits ${binfile}_dir/seg.fits ${boutfile}_dir/catalog.fits ${binfile}.psf_stars.fits psf_grid0.fits $npsf $x0 $y0 0
        ps1match.py ${binfile}_dir/catalog.fits ${boutfile}_dir/catalog.fits
        echo psf_deconvolve psf_grid.fits psf_grid0.fits 
        psf_deconvolve psf_grid.fits psf_grid0.fits 
        doutfile=${binfile}.diff.fits
        echo psf_subtract $infile ${binfile}.wt.fits $outfile ${boutfile}.wt.fits delta_psf_grid.fits \!${binfile}.diff.fits \!b${infile}.diff.wt.fits
        psf_subtract $infile ${binfile}.wt.fits $outfile ${boutfile}.wt.fits delta_psf_grid.fits \!$doutfile \!${binfile}.diff.wt.fits
        [ -d "${binfile}.diff_dir" ] && rm -r ${binfile}.diff_dir
        mkdir ${binfile}.diff_dir
        cp ${binfile}_dir/catalog.fits ${binfile}.diff_dir
        echo calc_psf $doutfile ${binfile}.diff.wt.fits ${binfile}_dir/seg.fits ${binfile}.diff_dir/catalog.fits ${binfile}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
        calc_psf $doutfile ${binfile}.diff.wt.fits ${binfile}_dir/seg.fits ${binfile}.diff_dir/catalog.fits ${binfile}.psf_stars.fits psf_grid.fits $npsf $x0 $y0 0
        cat2radec.py ${binfile}.diff_dir/catalog.fits ${binfile}.diff_dir/${binfile}.diff_radec.txt

        # now report detections
        mag0=`gethead MAG0 $infile`; vmag0=`gethead DMAG0 $infile | awk '{print $1*$1}'`
        [ "$mag0" ] || mag0=0.0
        [ "$vmag0" ] || vmag0=0.0
        cat0=${binfile}_dir/${binfile}_radec.txt
        cat1=${binfile}.diff_dir/${binfile}.diff_radec.txt
        cat2=${binfile}.diff_dir/${binfile}.diff_imsub.txt
        paste $cat1 $cat0 | awk '{dm=$3-$18; if(dm<0) dm=-dm; if ($3>0 && $NF<0 && $4<0.21715 && dm<sqrt(0.01+$4*$4)) printf("%d %.6f %.6f\n",$15,$3+'$mag0',sqrt($4*$4+'$vmag0'))}' | sort -nr -k 1 > $cat2

        [ -f "$doutfile" ] && fits2jpg.py $doutfile 1024 linvert
    fi

else
    echo "No files to stack."
fi

rm ps1_list$$.txt 2>/dev/null
