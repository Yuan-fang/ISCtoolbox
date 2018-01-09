#! /usr/bin/env bash

# Credit to Zarrar Shehzad
# https://github.com/czarrar/connectir/wiki/Generating-the-group-mask


# I/O Paths
funcpaths=($(cat func_paths.txt)) # inputs
maskdir=$(pwd)  # output directory for subject masks
mask_file="group_mask.nii.gz" # output file

# Create individual subject brain masks
n=${#funcpaths[@]}
for (( i = 0; i < $n; i++ )); do
    func=${funcpaths[$i]}
    mask="${maskdir}/mask${i}.nii.gz"
    fslmaths ${func} -Tstd -bin ${mask}
done

# Take the mean of the masks
# (i.e., proportion of subjects with value in each voxel)
3dMean -prefix group_prop_subjs.nii.gz ${maskdir}/mask*.nii.gz
  
# Note: you may also want to concatenate all the subject masks together to view them
# Here's how you could do that (optional!)
fslmerge -t ${maskdir}/subject_masks.nii.gz ${maskdir}/mask*.nii.gz

# Get voxels with all subjects having a value
3dcalc -a group_prop_subjs.nii.gz -expr 'equals(a,1)' -prefix ${mask_file}

# And if you want to mask the group mask further
grey_matter="/nfs/h1/workingshop/zhaoyuanfang/ISC-space/analysis/grpmask/greymatter.nii.gz"
3dcalc -a ${mask_file} -b ${grey_matter} -expr 'a*step(b)' -prefix grey_grpmask.nii.gz
