#!/bin/bash

file_list=$1
shift

bin=4
gain=6.2
maskfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 2 ; }
file0_wt=${file0%'.fits'}.wt.fits
[ -f $file0_wt ] || { echo "Cannot find first weight file $file0_wt" ; exit 3 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

cam=`gethead CCD_NAME $file0`

# resample the files in parallel (resampling preserves good seeing)
n1=`cat $file_list | wc -l | awk '{printf("%.0f\n",1+$1/'$NBATCH')}'`
rm sublist_$$_*.txt 2>/dev/null
awk 'BEGIN{i=0}{if(NR%'$n1'==0) i=i+1; print > "sublist_'$$'_"i".txt"}' $file_list

sargs3="-c ${SWARP_DIR}/coatli_redux.swarp -RESAMPLE Y -WEIGHT_SUFFIX .wt.fits -FSCALE_KEYWORD FSCALE0"

for file in `ls sublist_$$_*.txt`; do
    echo swarp @$file $sargs3 -COMBINE N -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits
    swarp @$file $sargs3 -COMBINE N -COMBINE_TYPE WEIGHTED -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits 2>/dev/null &
    gonogo
done
wait; go_iter=0

for file in `cat $file_list`; do
    base=${file%'.fits'}
    #[ -f ${base}.resamp.weight.fits ] && mv ${base}.resamp.weight.fits ${base}.resamp.wt.fits &
    fscale0=`gethead FSCALE0 ${base}.fits`
    sethead FSCALE0=$fscale0 ${base}.resamp.fits
    mv ${base}.resamp.fits ${base}.fits
    mv ${base}.resamp.weight.fits ${base}.wt.fits
done

n=`cat $file_list | wc -l`
m=`expr $n + 1`

file_list_iter=s_$$_$file_list
sort -n $file_list > $file_list_iter

awk '{printf("sethead IMID=%d IMNUM=1 IMPID=0 %s\n",NR,$1)}' $file_list_iter | sh

for file in `cat $file_list_iter`; do
    xmin=`gethead XMIN $file`
    [ "$xmin" ] || {
        nx=`gethead NAXIS1 $file`
        ny=`gethead NAXIS2 $file`
        sethead XMIN=1 XMAX=$nx YMIN=1 YMAX=$ny $file
    }
done

nh=`echo $n | awk '{print int($1/2)}'`
rm s?_$$_$file_list 2>/dev/null
if [ "$nh" -ge 1 ]; then
    sed -n "1,${nh}p" s_$$_$file_list > sa_$$_$file_list
    nh1=`expr $nh + 1`
    sed -n "${nh1},\$p" s_$$_$file_list > sb_$$_$file_list
fi

cp stack_${cam}.head stack_${cam}a.head
cp stack_${cam}.head stack_${cam}b.head

for file_list_iter in `ls s?_$$_$file_list`; do

    ab=`echo $file_list_iter | cut -c2`
    n=`cat $file_list_iter | wc -l`
    go=1

    while [ $go -eq 1 ]; do

        nmod=`echo $n $bin | awk '{print int($1/$2)/$1}'`
        awk '{i=int((NR-1)*'$nmod'); print > "nlist_'$$'_"i"_.txt"}' $file_list_iter
        n_next=`echo $n $nmod | awk '{print 1+int(($1-1)*$2)}'`

        ofile=
        if [ "$n_next" -eq 1 ]; then
            go=0
            ofile=stack_${cam}${ab}.fits
        fi

        for flist in `ls nlist_$$_*_.txt`; do
            sethead IMPID=$m @$flist
            restack_coatli.sh $flist maskfile=$maskfile imid=$m do_phot=no wcs_corr=no ofile=$ofile gain=$gain &
            m=`expr $m + 1`
            gonogo
        done > $file_list_iter
        wait; go_iter=0
        sort -n -k 1 $file_list_iter > flist$$.tmp
        mv flist$$.tmp $file_list_iter
        rm nlist_$$_*_.txt
        n=`cat $file_list_iter | wc -l`

    done

    sethead IMPID=0 stack_${cam}${ab}.fits

done

# make the final stack
ls stack_${cam}a.fits stack_${cam}b.fits > $file_list_iter
sethead IMPID=0 @$file_list_iter
restack_coatli.sh $file_list_iter imid=0 ofile=stack_${cam}.fits do_phot=no wcs_corr=no gain=$gain
sethead IMPID=-1 stack_${cam}.fits

rm s*_$$_$file_list sublist_$$_*.txt
