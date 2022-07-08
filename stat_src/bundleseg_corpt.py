import os
import sys
from glob import glob
import nibabel as nib
import numpy as np


#a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_est/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

a =  nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/true_'+sys.argv[1]+'.nii').get_fdata()

#b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()
b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_1/corpt_'+sys.argv[1]+'.nii').get_fdata()


#a = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_'+ sys.argv[1] + '.nii').get_fdata()
#b = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+ sys.argv[1] + '.nii').get_fdata()


err = b - a
pe =  ( err / a) * 100
ape = np.abs(pe)

print(np.nanmin(ape))
print(np.nanmax(ape))
print(err.shape)
#mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_mask.nii.gz'))
mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/tractsegout/TOM_trackings/*_mask.nii.gz'))
mask_files.sort()
values = []
for mask_file in mask_files:
    mask = nib.load(mask_file).get_fdata()
    roi = ape * mask
    mean_roi = np.nanmean(roi)
    #print(mask_file, mean_roi)
    values.append(mean_roi)

print(values)
print(min(values))
print(max(values))
