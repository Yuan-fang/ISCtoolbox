#! /usr/bin/env bash
base='/nfs/h1/workingshop/zhaoyuanfang/DPrest'
item='resting/002/sm6_inorm_bp0.01_0.1_csf_wm_gs_mc_confrm_lin_3mm.nii.gz'

for i in `cat sessid126`
do
    if [ ${i:0:1} == 'D' ]
    then
        echo ${base}/'restnii_DP'/${i}/${item} >> func_paths.txt
    else
        echo ${base}/'restnii_OldControl'/${i}/${item} >> func_paths.txt
    fi
done
