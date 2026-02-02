#!/bin/bash

file=`readlink -f $1`
shift

inst=$1
if [ "$inst" ]; then
   shift
else
    inst=coatli
fi
[ "$inst" ] || inst=coatli

do_dao=no
if [ "$1" = "dao" ]; then
    do_dao=yes
    shift
fi
do_force=no
if [ "$1" = "force" ]; then
    do_force=yes
    shift
fi

srcfile=`echo $1 | awk -F= '/srcfile/{print $2}'`
[ "$srcfile" ] && shift

extras=$@

base=${file%'.fits'}
rms_file=${base}.rms.fits
wt_file=${base}.wt.fits
[ -f "$wt_file" ] || wt_file=${base}.weight.fits

file0=`basename $file`
base=${file0%'.fits'}
base_dir=`dirname $file`
dir=${base}_dir
[ -d $dir ] || mkdir $dir

cd $dir

if [ "$srcfile" ]; then
    [ -f sky.list ] || radec2xy.py $file ${base_dir}/$srcfile > sky.list
fi
if [ -f sky.list ]; then
    echo "Running in assoc mode using sky.list file"
    inst=${inst}_assoc
fi

gain=`gethead GAIN $file`
[ "$gain" ] || gain=1.0

wfile=
if [ -f "$rms_file" ];then
    extras="-WEIGHT_TYPE MAP_RMS -WEIGHT_GAIN N -WEIGHT_IMAGE $rms_file -GAIN 1.e9 $extras"
    wfile=$rms_file
elif [ -f "$wt_file" ];then
    extras="-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_GAIN Y -WEIGHT_IMAGE $wt_file -WEIGHT_THRESH 0 -GAIN $gain $extras"
    wfile=$wt_file
else
    extras="-WEIGHT_TYPE NONE $extras"
fi

if [ "$do_dao" = "yes" ]; then
    wmap_file=${file%'.fits'}.wmap_mask.fits
    [ -f "$wmap_file" ] || extras="$extras -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME mask.fits"
fi

t0=`gethead DATE-OBS $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
t1=`gethead DATE-OBE $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
[ "$t1" ] || t1=$t0

ccd=`gethead CCD_NAME $file`
[ "$ccd" ] || ccd=`echo $file0 | cut -c16-17`

sat_level=`gethead SATURATE $file`
[ "$sat_level" ] || sat_level=60000.0
extras="-SATUR_LEVEL $sat_level $extras -CATALOG_NAME test.cat "

echo sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${inst}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw
sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${inst}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw

exposure=`gethead EXPTIME $file`
am=`gethead AIRMASS $file`
[ "$am" ] && airmass=`gethead AM_COEF $file | awk '{if('$am'>1) print $1*('$am'-1); else print 0.0}'`
[ "$airmass" ] || airmass=0.0

[ -f test_assoc.cat ] && mv test_assoc.cat test.cat

echo "# RA DEC mag dmag mag_big dmag_big FWHM x y xa ya x2a y2a expos num (ccd $ccd gain $gain exposure $exposure airmass_ext $airmass sex_zero 25.0 t0 $t0 t1 $t1 )" > ${base}_radec.txt

xmin=`gethead XMIN $file`
if [ "$xmin" ]; then
    xmax=`gethead XMAX $file`
    ymin=`gethead YMIN $file`
    ymax=`gethead YMAX $file`
    grep -v '#' test.cat | awk '{if($1>='$xmin'+5 && $1<='$xmax'-5 && $2>='$ymin'+5 && $2<='$ymax'-5) print}' > test.tmp
    mv test.tmp test.cat
fi

n=`grep -v '#' test.cat | wc -l`
if [ "$n" -gt 0 ]; then

    grep -v '#' test.cat | awk '{print $1,$2,$3}' > xy_in.txt
    grep -v '#' test.cat | awk '{if (NF==13) print $11,$12,$4,$6,$5,$7,$13,$1,$2; else print $25,$26,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$27,$1,$2}' >> ${base}_radec.txt

    # replace x and y positions with relative positions
    grep '#' ${base}_radec.txt > ${base}_radec.tmp
    grep -v '#' ${base}_radec.txt > ${base}_radec.tmp0
    fstr=`echo $file0 | cut -c1-2`
    if [ -f "$wt_file" -a "$fstr" != "20" ]; then
        imxy2expos.py $wt_file xy_in.txt > xy_out.txt 2>/dev/null
    else
        sx=`gethead NAXIS1 $file`; sy=`gethead NAXIS2 $file`
        awk '{x=$1/'$sx';y=$2/'$sy'; print x,y,0.5*(3*x*x-1),0.5*(3*y*y-1),'$exposure',$3}' xy_in.txt > xy_out.txt
    fi
    paste -d ' ' ${base}_radec.tmp0 xy_out.txt >> ${base}_radec.tmp
    grep -v 'nan' ${base}_radec.tmp > ${base}_radec.txt
    rm ${base}_radec.tmp0 ${base}_radec.tmp

    if [ "$do_dao" = "yes" ]; then
        cp ${base}_radec.txt ${base}_radec0.txt
        if [ -f "${base_dir}/force_phot.txt" -a "$do_force" = "yes" ]; then
	    radec2xy.py $file ${base_dir}/force_phot.txt | awk '{print $1,$2,-NR}' > xy_in.txt
            if [ -s xy_in.txt ]; then
                echo "Forcing photometry at positions in ${base_dir}/force_phot.txt:"
                cat xy_in.txt
                if [ -f "$wt_file" -a "$fstr" != "20" ]; then
                    imxy2expos.py $wt_file xy_in.txt > xy_out.txt 2>/dev/null
                else
                    sx=`gethead NAXIS1 $file`; sy=`gethead NAXIS2 $file`
                    awk '{x=$1/'$sx';y=$2/'$sy'; print x,y,0.5*(3*x*x-1),0.5*(3*y*y-1),'$exposure',$3}' xy_in.txt > xy_out.txt
                fi
                dat=`grep -v '#' ${base}_radec.txt | awk '{printf("%f,%f,%f,%f,%f\n",$3,$4,$5,$6,$7);exit}'`
                [ "$dat" ] || dat="19.0,0.1,19.0,0.1,10.0"
		awk '{print $1,$2}' xy_in.txt > xy_in0.txt
                paste ${base_dir}/force_phot.txt xy_in0.txt xy_out.txt | eval awk "'{print \$1,\$2,$dat,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}'" >> ${base}_radec.txt
		rm xy_in0.txt
            fi
        fi
        fwhm=`gethead FWHM $file`
        [ "$fwhm" ] || {
            fwhm=`quick_mode test.cat n=NF min=0`
            sethead FWHM=$fwhm $file
        }
        [ -f "$wmap_file" ] || {
            echo sources2mask.py $wt_file mask.fits $wmap_file $fwhm
            sources2mask.py $wt_file mask.fits $wmap_file $fwhm
        }
        echo psf_fit.py $file ${base}_radec.txt $psfonly --satlevel $sat_level $aperture
        psf_fit.py $file ${base}_radec.txt $psfonly --satlevel $sat_level $aperture > ${base}_radec1.txt
        [ -f ${base}_psf_stats.txt ] && mv ${base}_psf_stats.txt psf_stats.txt
        if [ -s ${base}_radec1.txt ]; then
            mv ${base}_radec1.txt ${base}_radec.txt
        else
            echo "PSF fitting not applied or failed."
        fi
        fwhm_eff=`grep FWHMeff ${base}_radec.txt | awk '{print $NF}'`
        ap_corr=`grep FWHMeff ${base}_radec.txt | awk '{print $(NF-2)}'`
        echo "FWHMeff: $fwhm_eff"
        [ "$fwhm_eff" ] && sethead FWHMeff=$fwhm_eff APCorr=$ap_corr $file
    fi

fi
