#!/bin/bash

infile=$1
shift

catdir=.
outfile=ps1.fits
maskfile=mask.fits
setup_only=no
do_wait=no

local_stack_dir=$REDUX_BASE_DIR/ps1_stacks

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f "$infile" ] || { echo "Input file $infile does not exist." ; exit 1 ; }

ps1_files=`get_ps1imagelist.py 0 0 1 f $infile`
echo "ps1_files: $ps1_files"

[ "$ps1_files" ] || { echo "Cannot find ps1 files to download." ; exit 2 ; }

filter=`gethead FILTER $infile`
[ "$filter" = "B" ] && filter=g
[ "$filter" = "gri" ] && filter=r
[ "$filter" = "zy" ] && filter=z
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
    lfile=$local_stack_dir/$cell/$scell/${file}.fz
    if [ -f "$ofile" ]; then
        echo "$ofile exists from previous run."
    else
        [ -f $lfile ] && cp $lfile $ofile
    fi
    if [ -f "$ofile" ]; then
        mkdir -p $local_stack_dir/$cell/$scell 2>/dev/null
        #[ -f $lfile ] || cp $ofile $lfile
    else
        echo "$ofile" >> ps1_downloading$$.txt
        wget "${url1}/$file" $wg -O $ofile >/dev/null 2>&1 &
    fi

    base=${file%'.fits'}
    ofile=${pdir}/${base}.var.fits.fz
    lfile=$local_stack_dir/$cell/$scell/${base}.var.fits.fz
    if [ -f "$ofile" ]; then
        echo "$ofile exists from previous run."
    else
        [ -f $lfile ] && cp $lfile $ofile
    fi
    if [ -f "$ofile" ]; then
        mkdir -p $local_stack_dir/$cell/$scell 2>/dev/null
        #[ -f $lfile ] || cp $ofile $lfile
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

    sargs="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE Y -SUBTRACT_BACK N -WEIGHT_SUFFIX .wt.fits -WEIGHT_TYPE MAP_WEIGHT -BLANK_BADPIXELS N"
    echo swarp @ps1_list$$.txt $sargs -COMBINE_TYPE AVERAGE -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME $outwfile
    swarp @ps1_list$$.txt $sargs -COMBINE_TYPE AVERAGE -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME $outwfile
    rm ${boutfile}.head

    gain=`gethead GAIN $outfile`
    sethead GAIN=1.E9 $outfile  # this file has an rms weightmap

    # background subtraction from sextractor
    run_sex.sh $outfile coatli -DETECT_THRESH 20.0 -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME back.fits
    fwhm=`quick_mode ${boutfile}_dir/test.cat n=NF min=0`
    [ "$fwhm" ] || fwhm=2.0
    phot_diam=`echo $fwhm | awk '{printf("%.1f\n",1+int($1*3))}'`
    sethead FWHM=$fwhm PHOTDIAM=$phot_diam $outfile
    imarith $outfile ${boutfile}_dir/back.fits sub tmp/$outfile
    mv tmp/$outfile $outfile

    cp $outfile n$outfile       # we will create a standard (non-rms) weightmap for this
    sethead GAIN=$gain n$outfile

    # make a weightmap appropriate for psf fitting
    if [ -f "$maskfile" ]; then
        ps1wt2wt.py $outwfile $maskfile n$outwfile
        var0=`gethead VAR0 n$outwfile`
        sethead VAR0=$var0 n$outfile
    fi

else
    echo "No files to stack."
fi

rm ps1_list$$.txt 2>/dev/null
