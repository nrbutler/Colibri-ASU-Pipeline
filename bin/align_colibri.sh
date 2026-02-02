#!/bin/bash

file_list=$1
shift

nstars_min=5
inst=coatli
thresh=10.0

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

echo find_stars.sh $file_list inst=$inst nstars_min=$nstars_min thresh=$thresh
find_stars.sh $file_list inst=$inst nstars_min=$nstars_min thresh=$thresh
mv stars_$file_list $file_list
nfiles=`cat $file_list | wc -l`

align_ref_file=align_file.txt
if [ "$nfiles" -gt 0 ]; then

    # make sure an existing align reference file is usable
    if [ -s "$align_ref_file" ]; then
        afile=`head -1 $align_ref_file | awk '{print $2}'`
        dt0=`gethead EXPTIME $afile | awk '{printf("%.0f\n",$1)}'`
        dtmax=`gethead -a EXPTIME @$file_list | sort -nr -k 2 | awk '{printf("%.0f\n",$2);exit}'`
        if [ "$dt0" -lt "$dtmax" ]; then
           echo "Removing alignment reference file $align_ref_file ($dt0<$dtmax)"
           rm $align_ref_file
           sethead MATCHED=0 @$file_list
        fi
    fi

    if [ -s "$align_ref_file" ]; then
        echo "Keeping alignment reference file from a previous run."
    else
        #afile=`gethead -a NSTARS @$file_list | sort -n -k 2 | awk '{if(NR==int('$nfiles'/2)+1) print $1}'`
        afile=`gethead -a NSTARS @$file_list | sort -nr -k 2 | awk '{if(NR>'$nfiles'/10) {print $1; exit}}'`
        sethead MATCHED=1 $afile
        echo "# $afile" > $align_ref_file
        if [ -f "$afile" ]; then
            base=${afile%'.fits'}
            cat ${base}_dir/${base}_radec.txt >> $align_ref_file
        fi
    fi
fi

if [ -s "$align_ref_file" ]; then

    rm report$$.txt 2>/dev/null
    afile=`head -1 $align_ref_file | awk '{print $2}'`
    echo "Aligning images to $afile"
    nfiles1=`cat $file_list | grep -v $afile | wc -l`
    for file in `cat $file_list | grep -v $afile`; do
        base=${file%'.fits'}
        sexlist_cc.py $file ${base}_dir/${base}_radec.txt $align_ref_file X 0 0 0 25 $nstars_min >> report$$.txt 2>/dev/null &
        gonogo
    done
    wait; go_iter=0

    nmatched0=`gethead NMATCHED $afile`
    if [ "$nmatched0" -gt 0 ]; then
        echo "Alignment reference has about $nmatched0 matchable stars"
    else
        nmatched0=`gethead -a NMATCHED @$file_list | grep -v '___' | sort -nr -k 2 | awk '{print $2;exit}'`
        sethead NMATCHED=$nmatched0 FSCALE0=1.0 $afile
    fi

    nmatch=`cat report$$.txt 2>/dev/null | wc -l`
    echo "Matched $nmatch of $nfiles1 files..."
    cat report$$.txt 2>/dev/null

    nmatch_med=`gethead -a NMATCHED @$file_list | sort -n -k 2 | awk '{if(NR==int('$nfiles'/2)+1) print $2}'`

    grep $afile $file_list > ${file_list}.tmp
    awk '{print $5,$6}' report$$.txt 2>/dev/null | awk '{if($1>'$nmatch_med'/4) print $2}' | awk -F/ '{print $2}' | sed -e 's/_radec\.txt/\.fits/g' >> ${file_list}.tmp
    #awk '{print $6}' report$$.txt | awk -F/ '{print $2}' | sed -e 's/_radec\.txt/\.fits/g' >> ${file_list}.tmp
    mv ${file_list}.tmp $file_list
    sethead MATCHED=1 @$file_list
    rm report$$.txt 2>/dev/null

fi
