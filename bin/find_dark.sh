#!/bin/sh

day=$1
shift
ccd=$1
shift
bin=$1
shift
rmode=$1
shift

dir=${REDUX_BASE_DIR}/dark_bank

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -d "$dir" ] || mkdir $dir

cd $dir

ls dark*.fits* 2>/dev/null > file_list0$$.txt
sed -e 's/_/ /g' -e 's/\.fits/ /g' file_list0$$.txt | awk '{print $3,$4,$5,$6}' > file_list1$$.txt
paste file_list0$$.txt file_list1$$.txt | eval "awk '{if(\$2~/$ccd/ && \$4==$bin && \$5~/$rmode/) print \$1}'" > file_list$$.txt
rm file_list0$$.txt file_list1$$.txt

cut -c6-13 file_list$$.txt > days$$.txt
paste days$$.txt file_list$$.txt > days$$.tmp
sort -n -k 1 days$$.tmp > days$$.txt
rm days$$.tmp

day0=`awk '{if($1>'$day') exit; day=$1}END{print day}' days$$.txt`
[ "$day0" ] || day0=`awk '{print $1; exit}' days$$.txt`

for file in `awk '{if($1~/'$day0'/) print $2}' days$$.txt`; do
    echo ${dir}/$file
done
rm days$$.txt file_list$$.txt
