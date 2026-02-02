#!/bin/bash

file_list=$1
shift

nstars_min=10
fwhm_fac=3  #  exclude images with fwhm > fwhm_fac * median(fwhm)
thresh=1.5
cal_mag_max=18.0
inst=coatli
do_dao=no
sex_args=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

[ "$thresh" ] || thresh=1.5

dao=
if [ "$do_dao" = "yes" ]; then
    dao=dao
fi

for file in `cat $file_list`; do
    base=${file%'.fits'}
    if [ -d "${base}_dir" ]; then
        echo "${base}_dir exists, not remaking" > /dev/null 2>&1
    else
        run_sex.sh $file $inst $dao -DETECT_THRESH $thresh $sex_args > /dev/null 2>&1 &
        gonogo
    fi
done
wait; go_iter=0

function get_file_info() {
    local file=$1
    local base=${file%'.fits'}
    local nstars=`grep -v '#' ${base}_dir/${base}_radec.txt | wc -l`
    if [ "$nstars" -gt 0 ]; then
        local fwhm=`quick_mode ${base}_dir/${base}_radec.txt n=7`
        if [ "$fwhm" ]; then
            head -1 ${base}_dir/${base}_radec.txt > ${base}_dir/${base}_radec.tmp
            grep -v '#' ${base}_dir/${base}_radec.txt | awk '{if($7>'$fwhm'/2. && $7<'$fwhm'*4) print}' >> ${base}_dir/${base}_radec.tmp
            mv ${base}_dir/${base}_radec.tmp ${base}_dir/${base}_radec.txt
            nstars=`grep -v '#' ${base}_dir/${base}_radec.txt | wc -l`
            echo $file $fwhm $nstars
        fi
    fi
    [ "$fwhm" ] || fwhm=0.0
    sethead NSTARS=$nstars FWHM=$fwhm $file
}

for file in `cat $file_list`; do
    get_file_info $file &
    gonogo
done > fwhm$$.txt
wait; go_iter=0

fwhm0=`quick_mode fwhm$$.txt n=2`
if [ "$fwhm0" ]; then
    echo "Median FWHM of stars in images is $fwhm0 pixels"
    awk '{if($3>'$nstars_min' && $2<'$fwhm_fac'*'$fwhm0') print $1}' fwhm$$.txt > stars_$file_list
else
    awk '{if($3>'$nstars_min') print $1}' fwhm$$.txt > stars_$file_list
fi
rm fwhm$$.txt

n0=`cat $file_list | wc -l`
n1=`cat stars_$file_list | wc -l`
echo "Keeping $n1 of $n0 images with $nstars_min or more stars."
