#!/bin/bash

bias=$1
shift
infile_list=$1
shift

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ "$NBATCH" ] || NBATCH=`cat /proc/cpuinfo | grep processor | wc -l`

[ -f $bias ] || { echo "no file $bias" ; exit 3 ; }

[ -d tmp ] || mkdir tmp

echo "Subtracting frame $bias from:"
cat $infile_list

function doit() {
    local file=$1
    imarith $file $bias sub tmp/$file
    mv tmp/$file $file
}

n=0
for file in `cat $infile_list`; do
    doit $file &
    n=`expr $n + 1`
    go=`echo $n | awk '{print $1 % '$NBATCH'}'`
    [ "$go" -eq 0 ] && wait
done
wait
