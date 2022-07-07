import os
import sys
from glob import glob
import nibabel as nib
import numpy as np


a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_est/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

err = b - a
pe =  ( err / a) * 100
ape = np.abs(pe)

mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_mask.nii.gz'))

for mask_file in mask_files:
    mask = nib.load(mask_file).get_fdata()
    roi = ape * mask
    mean_roi = np.nanmean(roi)
    print(mask_file, mean_roi)
