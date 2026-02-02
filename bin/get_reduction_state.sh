#!/bin/bash

#cd $REDUX_WORK_BASE_DIR
cd /home/colibri

base=_usr_local_var_colibri_
for file in `ls ${base}*/redux_colibri_C*/current_state.txt`; do
    state=`cat $file`
    obs=`dirname $file | tr '_' ' ' | awk '{print $5,$6,$7,$8,$12,$13}'`
    echo $state $obs
done
