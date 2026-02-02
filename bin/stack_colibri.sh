#!/bin/bash

file_list=$1
shift

#matching,wcs
mask_sigma=6.0
bg_order=4

catdir=

# base reduction
ps=0.000195
inst=coatli
gain=6.2

do_nir_sky=yes
new_mask=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 2 ; }
file1=`sort $file_list | tail -1`

[ "$catdir" ] || { echo "No catalog dir specified" ; exit 3 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

cam=`basename $file0 | cut -c16-17`
nx=`gethead NAXIS1 $file0`
ny=`gethead NAXIS2 $file0`
cx=`gethead CRPIX1 $file0`
cy=`gethead CRPIX2 $file0`

nstack=`cat $file_list | wc -l`

if [ "$nstack" -gt 0 ]; then

    # make a new mask if pointing has shifted or depth has increased markedly
    # exptime current > 0.8 exp_time_last
    # pointing not changed (like never true until end)
    if [ -f "stack_${cam}.head" -a "$new_mask" = "no" ]; then
        echo "Using exisining stack_${cam}.head for stack organization."
    else
        sargs0="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLING_TYPE NEAREST -SUBTRACT_BACK N -COMBINE N -HEADER_ONLY Y"
        sargs0="$sargs0 -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits -WEIGHT_TYPE NONE"
        rm stack_${cam}.head 2>/dev/null
        echo swarp @$file_list $sargs0 -MEM_MAX 1024 -COMBINE_BUFSIZE 1024
        swarp @$file_list $sargs0 -MEM_MAX 1024 -COMBINE_BUFSIZE 1024
        sethead CD1_1=-$ps CD2_2=$ps stack_${cam}.fits

        crpix_x=`gethead CRPIX1 stack_${cam}.fits | awk '{printf("%.0f\n",$1)}'`
        crpix_x=`echo $crpix_x | awk '{printf("%.1f\n",$1-0.5)}'`
        crpix_y=`gethead CRPIX2 stack_${cam}.fits | awk '{printf("%.0f\n",$1)}'`
        crpix_y=`echo $crpix_y | awk '{printf("%.1f\n",$1-0.5)}'`
        sethead CRPIX1=$crpix_x CRPIX2=$crpix_y stack_${cam}.fits

        imhead -f stack_${cam}.fits > stack_${cam}.head
        sethead BITPIX=0 stack_${cam}.head
    fi

    # now we align the new matched files
    ra0=`gethead CRVAL1 stack_${cam}.head`
    dec0=`gethead CRVAL2 stack_${cam}.head`
    function copy_files() {
        local file=$1
        cp $file f$file
        local base=${file%'.fits'}
        local wfile=${base}.wt.fits
        cp $wfile f$wfile
        imhead -f $file > f${base}.head1
        local ra1=`gethead CRVAL1 $file`; dec1=`gethead CRVAL2 $file`
        local dx=`echo $ra1 $ra0 | awk '{printf("%.0f\n",($2-$1)/'$ps')}'`
        local dy=`echo $dec1 $dec0 | awk '{printf("%.0f\n",($2-$1)/'$ps')}'`
        local ra2=`echo $ra1 $dx | awk '{printf("%.8f\n",$1+$2*'$ps')}'`
        local dec2=`echo $dec1 $dy | awk '{printf("%.8f\n",$1+$2*'$ps')}'`
        local cx1=`echo $cx $dx | awk '{print $1-$2}'`
        local cy1=`echo $cy $dy | awk '{print $1+$2}'`
        local sky=`gethead SKYLEV $file`
        sethead CRVAL1=$ra2 CRVAL2=$dec2 CRPIX1=$cx1 CRPIX2=$cy1 f$file
        sethead BZERO=-$sky f$file
        echo f$file
    }

    for file in `cat $file_list`; do
        copy_files $file &
        gonogo
    done > f$file_list
    wait; go_iter=0

    sort -n f$file_list > f${file_list}.tmp
    mv f${file_list}.tmp f$file_list

    # record some stats
    sat_level=`gethead -a SATURATE @f$file_list | awk 'BEGIN{s=1.e56}{if($2<s) s=$2}END{print s}'`
    exptime=`gethead -a EXPTIME @f$file_list | awk '{s=s+$2}END{printf("%f\n", s)}'`
    #igain=`echo $gain $exptime | awk '{print $1*$2}'`    # weighted stack gain
    wgain=`gethead -a GAIN VAR0 @f$file_list | awk '{s=s+1/($2*$3*$3);s0=s0+1/$3}END{print s0*s0/s}'` # better version
    var0=`gethead -a VAR0 @f$file_list | awk '{s=s+1/$2}END{printf("%.6f\n",1/s)}'`
    sky0=`gethead -a SKYLEV @f$file_list | awk '{n=n+1;w0=w0+1/$2}END{s=n/w0; if(s<0) s=0; print s}'`

    # we need to separate files by exposure time and mount rotation
    gethead EXPTIME SMTMRO @f$file_list | awk '{print $1>"f'$file_list'_"$2"-"$3".txt"}'
    nmax=`ls f${file_list}_*.txt | wc -l`
    ls f${file_list}_*.txt

    maskfile=maskstack_${cam}.fits
    if [ -f "$maskfile" -a "$new_mask" = "no" ]; then
        echo "Using exising mask $maskfile"
    else
        maskfile=maskstack0_${cam}.fits

        # build a stack for source masking
        sargs1="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE N -RESAMPLE_SUFFIX .fits -SUBTRACT_BACK Y -WEIGHT_SUFFIX .wt.fits"
        if [ "$nmax" -gt 1 ]; then
            n=0
            rm list$$.txt 2>/dev/null
            for flist in `ls f${file_list}_*.txt`; do
                echo stack0_${cam}_${n}.fits >> list$$.txt
                echo swarp @$flist $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack0_${cam}_${n}.fits -WEIGHTOUT_NAME stack0_${cam}_${n}.wt.fits
                swarp @$flist $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack0_${cam}_${n}.fits -WEIGHTOUT_NAME stack0_${cam}_${n}.wt.fits 2>/dev/null &
                gonogo
                n=`expr $n + 1`
            done
            wait; go_iter=0
        else
            cp f${file_list}_*.txt list$$.txt
        fi
        cp stack_${cam}.head stack0_${cam}.head
        echo swarp @list$$.txt $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack0_${cam}.fits -WEIGHTOUT_NAME stack0_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024
        swarp @list$$.txt $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack0_${cam}.fits -WEIGHTOUT_NAME stack0_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024 2>/dev/null
        rm stack0_${cam}.head list$$.txt 2>/dev/null
        sethead SKYLEV=$sky0 SATURATE=$sat_level EXPTIME=$exptime VAR0=$var0 GAIN=$wgain stack0_${cam}.fits

        # make a source mask
        run_sex.sh stack0_${cam}.fits $inst -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME mask.fits -BACK_SIZE 8
        echo sources2mask.py stack0_${cam}.wt.fits stack0_${cam}_dir/mask.fits $maskfile $mask_sigma
        sources2mask.py stack0_${cam}.wt.fits stack0_${cam}_dir/mask.fits $maskfile $mask_sigma
    fi
    echo $maskfile > maskfile_used.txt

    # register that mask to every file
    x0=`gethead CRPIX1 stack0_${cam}.fits`; y0=`gethead CRPIX2 stack0_${cam}.fits`
    nx0=`gethead NAXIS1 stack0_${cam}.fits`; ny0=`gethead NAXIS2 stack0_${cam}.fits`
    function sub_back() {
        local file=$1
        local base=${file%'.fits'}
        local xr=`gethead CRPIX1 $file | awk '{x=int($1-('$x0'));if ('$nx'-x>'$nx0') x='$nx'-'$nx0';printf("%.0f-%.0f\n",1-x,'$nx'-x)}'`
        local yr=`gethead CRPIX2 $file | awk '{y=int($1-('$y0'));if ('$ny'-y>'$ny0') y='$ny'-'$ny0';printf("%.0f-%.0f\n",1-y,'$ny'-y)}'`
        echo getfits $maskfile $xr $yr -o ${base}.wmap_mask.fits
        getfits $maskfile $xr $yr -o ${base}.wmap_mask.fits
        backsub.py $file ${base}.wmap_mask.fits $bg_order
    }

    for file in `cat f$file_list`; do
        sub_back $file &
        gonogo
    done
    wait; go_iter=0

    # determine an on-chip background
    rm good_lfiles.txt bad_lfiles.txt 2>/dev/null
    for lfile in `ls f${file_list}_*.txt`; do
        gethead -a CRPIX1 CRPIX2 @$lfile > xy$$.txt
        x0=`awk '{s=s+$2}END{print s/NR}' xy$$.txt`
        y0=`awk '{s=s+$3}END{print s/NR}' xy$$.txt`
        vx=`awk '{s=s+($2-('$x0'))*($2-('$x0'))}END{print s/NR}' xy$$.txt`
        vy=`awk '{s=s+($3-('$y0'))*($3-('$y0'))}END{print s/NR}' xy$$.txt`
        rms=`echo $vx $vy | awk '{printf("%.0f\n",sqrt($1+$2))}'`
        if [ "$rms" -gt 10 ]; then
            echo $lfile >> good_lfiles.txt
        else
            echo $lfile >> bad_lfiles.txt
        fi
    done
    rm xy$$.txt 2>/dev/null

    [ "$do_nir_sky" = "no" ] && rm good_lfiles.txt 2>/dev/null

    # scale by the sky level
    for lfile in `cat good_lfiles.txt 2>/dev/null`; do
        gethead -a SKYLEV @$lfile | awk '{s=$2; if(s<1) s=1.; s0='$sky0'; if(s0<1) s0=1; printf("sethead BSCALE=%.6f %s\n",s0/s,$1)}' | sh
    done

    sargs2="-c ${SWARP_DIR}/coatli_redux.swarp -COMBINE_TYPE MEDIAN -RESAMPLE N -HEADER_SUFFIX .head1 -RESAMPLE_SUFFIX .fits -SUBTRACT_BACK N -WEIGHT_SUFFIX .wmap_mask.fits"
    n=0
    for lfile in `cat good_lfiles.txt 2>/dev/null`; do
        echo swarp @$lfile $sargs2 -IMAGEOUT_NAME back_${n}.fits -WEIGHTOUT_NAME back_${n}.wt.fits
        swarp @$lfile $sargs2 -IMAGEOUT_NAME back_${n}.fits -WEIGHTOUT_NAME back_${n}.wt.fits 2>/dev/null &
        gonogo
        n=`expr $n + 1`
    done
    wait; go_iter=0

    n=0
    for lfile in `cat good_lfiles.txt 2>/dev/null`; do
        subtract_images.sh back_${n}.fits $lfile &
        gonogo
        n=`expr $n + 1`
    done
    wait; go_iter=0

    for lfile in `cat bad_lfiles.txt 2>/dev/null`; do
        for file in `cat $lfile`; do
            backsub.py $file ${file%'.fits'}.wmap_mask.fits -1 &
            gonogo
        done
    done
    wait; go_iter=0

    # unscale by the sky level
    for lfile in `cat good_lfiles.txt 2>/dev/null`; do
        sethead BSCALE=1.0 @$lfile
    done

    # now create a median stack that will form the basis for a clipped stack
    if [ "$nmax" -gt 1 ]; then
        n=0
        rm list$$.txt 2>/dev/null
        for flist in `ls f${file_list}_*.txt`; do
            echo med_stack_${cam}_${n}.fits >> list$$.txt
            echo swarp @$flist $sargs1 -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME med_stack_${cam}_${n}.fits -WEIGHTOUT_NAME med_stack_${cam}_${n}.wt.fits
            swarp @$flist $sargs1 -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME med_stack_${cam}_${n}.fits -WEIGHTOUT_NAME med_stack_${cam}_${n}.wt.fits 2>/dev/null &
            gonogo
            n=`expr $n + 1`
        done
        wait; go_iter=0
    else
        cp f${file_list}_*.txt list$$.txt
    fi
    cp stack_${cam}.head med_stack_${cam}.head
    if [ "$nmax" -gt 1 ]; then
        echo swarp @list$$.txt $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME med_stack_${cam}.fits -WEIGHTOUT_NAME med_stack_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024
        swarp @list$$.txt $sargs1 -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME med_stack_${cam}.fits -WEIGHTOUT_NAME med_stack_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024 2>/dev/null
    else
        echo swarp @list$$.txt $sargs1 -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME med_stack_${cam}.fits -WEIGHTOUT_NAME med_stack_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024
        swarp @list$$.txt $sargs1 -COMBINE_TYPE MEDIAN -IMAGEOUT_NAME med_stack_${cam}.fits -WEIGHTOUT_NAME med_stack_${cam}.wt.fits -MEM_MAX 1024 -COMBINE_BUFSIZE 1024 2>/dev/null
    fi
    rm med_stack_${cam}.head list$$.txt 2>/dev/null

    sethead SKYLEV=$sky0 SATURATE=$sat_level EXPTIME=$exptime VAR0=$var0 GAIN=$wgain med_stack_${cam}.fits

    x0=`gethead CRPIX1 med_stack_${cam}.fits`; y0=`gethead CRPIX2 med_stack_${cam}.fits`
    nx0=`gethead NAXIS1 med_stack_${cam}.fits`; ny0=`gethead NAXIS2 med_stack_${cam}.fits`
    function do_weights() {
        local file=$1
        local fscale0=$2
        local base=${file%'.fits'}
        local wfile=${base}.wt.fits
        local wmfile=${base}.wmap_mask.fits
        local r0=`gethead CRVAL1 $file`
        local d0=`gethead CRVAL2 $file`
        local x=`gethead CRPIX1 $file`
        local y=`gethead CRPIX2 $file`
        local dx=`echo $x | awk '{printf("%.0f\n",'$x0'-$1+1)}'`
        local dy=`echo $y | awk '{printf("%.0f\n",'$y0'-$1+1)}'`
        local x1=`gethead NAXIS1 $file | awk '{printf("%.0f\n",'$dx'+$1-1)}'`
        local y1=`gethead NAXIS2 $file | awk '{printf("%.0f\n",'$dy'+$1-1)}'`
        [ "$x1" -gt "$nx0" ] && {
            dx=`echo $dx | awk '{printf("%.0f\n",$1-('$x1')+'$nx0')}'`
            x1=$nx0
        }
        [ "$y1" -gt "$ny0" ] && {
            dy=`echo $dy | awk '{printf("%.0f\n",$1-('$y1')+'$ny0')}'`
            y1=$ny0
        }
        echo weight_clip $file $wfile $wmfile med_stack_${cam}.fits[$dx:$x1,$dy:$y1] med_stack_${cam}.wt.fits[$dx:$x1,$dy:$y1] $fscale0
        weight_clip $file $wfile $wmfile med_stack_${cam}.fits[$dx:$x1,$dy:$y1] med_stack_${cam}.wt.fits[$dx:$x1,$dy:$y1] $fscale0
        sethead CRPIX1=$x CRPIX2=$y CRVAL1=$r0 CRVAL2=$d0 $file
    }

    nfiles=`cat f$file_list | wc -l`
    fscale0=`gethead -a FSCALE0 @f$file_list | sort -n -k 2 | awk '{if(NR==int('$nfiles'/2)) print $2}'`
    [ "$fscale0" ] || fscale0=1.0
    for file in `cat f$file_list`; do
        do_weights $file $fscale0 &
        gonogo
    done
    wait; go_iter=0

    # now produce all the weighted stacks
    rm -r stack_20* 2>/dev/null
    iterstack_coatli.sh f$file_list gain=$gain maskfile=$maskfile
fi

#
# finally set the wcs information
#
[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
[ -f stack_${cam}.rms.fits ] && rm stack_${cam}.rms.fits
run_astnet_blind.sh stack_${cam}.fits thresh=20.0
fwhm0=`quick_mode stack_${cam}_dir/stack_${cam}_radec.txt n=7 min=0`
[ "$fwhm0" ] || fwhm0=2.0
phot_diam=`echo $fwhm0 | awk '{printf("%.1f\n",1+int($1*3))}'`
sethead PHOTDIAM=$phot_diam FWHM=$fwhm0 stack_${cam}.fits

#
# possibly tune up the wcs
#
here=`pwd`
tunefile=gaia_radec.fits
if [ -s "$tunefile" ]; then
    echo "WCS tuning file $tunefile exists."
else
    if [ -f "${catdir}/gaia_radec_${cam}.fits" ]; then
        echo "Using ${catdir}/gaia_radec_${cam}.fits"
    else
        filef=`readlink -f stack_${cam}.fits`
        cd $catdir
        echo grab_gaia_local.sh $filef
        grab_gaia_local.sh $filef
        cd $here
    fi
    ln -s ${catdir}/gaia_radec_${cam}.fits gaia_radec.fits 2>/dev/null
fi

if [ -s "$tunefile" ]; then
   echo tune_wcs.py stack_${cam}.fits stack_${cam}_dir/stack_${cam}_radec.txt $tunefile
   tune_wcs.py stack_${cam}.fits stack_${cam}_dir/stack_${cam}_radec.txt $tunefile
fi

# try to get a matching ps1 stack
# get and process a ps1 stack for image differencing
filter=`gethead FILTER stack_${cam}.fits`
pfilter=$filter
[ "$pfilter" = "B" ] && pfilter=g
[ "$pfilter" = "gri" ] && pfilter=r
[ "$pfilter" = "zy" ] && pfilter=z
ps1_file=ps1_stack_${pfilter}.fits
echo get_ps1_ddoti.sh stack_${cam}.fits outfile=$ps1_file catdir=$catdir maskfile=$maskfile
get_ps1_ddoti.sh stack_${cam}.fits outfile=$ps1_file catdir=$catdir maskfile=$maskfile

# make sure trigger info is in the stack
pos_err0=`gethead ALUN $file1 | awk '{printf("%.6f\n",$1)}'`
if [ "$pos_err0" ]; then
    ra0=`gethead ALRA $file1 | awk '{printf("%.6f\n",$1)}'`
    dec0=`gethead ALDE $file1 | awk '{printf("%.6f\n",$1)}'`
    trig_time=`gethead ALEVT $file1`
    sethead ALUN=$pos_err0 ALRA=$ra0 ALDE=$dec0 ALEVT=$trig_time stack_${cam}.fits
fi
