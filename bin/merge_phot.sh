#!/bin/sh

bands=`ls -d photometry_20* | awk -F_ '{print $(NF-1)}' | sort -u`
nbands=`echo $bands | wc -w`

if [ "$nbands" -gt 1 ]; then

    for band in $bands ; do
        pdir=`ls -td photometry_20*_${band}_* | head -1`
        echo $pdir $band
        grep -v '#' ${pdir}/stack_20*_${band}_radec.txt
    done

fi
