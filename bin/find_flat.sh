#!/bin/sh

day=$1
shift
ccd=$1
shift
filter=$1
shift
bin=$1
shift

dir=${REDUX_BASE_DIR}/flat_bank

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -d "$dir" ] || mkdir $dir

cd $dir

ls flat*.fits* 2>/dev/null > file_list0$$.txt
sed -e 's/_/ /g' -e 's/\.fits/ /g' file_list0$$.txt | awk '{print $3,$4,$5}' > file_list1$$.txt
paste file_list0$$.txt file_list1$$.txt | eval "awk '{if(\$2~/$ccd/ && \$3~/$filter/ && \$4==$bin) print \$1}'" > file_list$$.txt
rm file_list0$$.txt file_list1$$.txt

if [ -s "file_list$$.txt" ]; then

    cut -c6-13 file_list$$.txt > days$$.txt
    paste days$$.txt file_list$$.txt > days$$.tmp
    sort -n -k 1 days$$.tmp > days$$.txt
    rm days$$.tmp

    file=`awk '{if($1>'$day') exit; file=$2}END{print file}' days$$.txt`
    [ "$file" ] || file=`awk '{print $2; exit}' days$$.txt`
    echo ${dir}/$file
    rm days$$.txt

fi
rm file_list$$.txt 2>/dev/null
