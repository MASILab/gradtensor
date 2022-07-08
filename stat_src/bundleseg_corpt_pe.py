import os
import sys
from glob import glob
import nibabel as nib
import numpy as np


#a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_est/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

#b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

a =  nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/true_'+sys.argv[1]+'.nii').get_fdata()

b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_1/corpt_'+sys.argv[1]+'.nii').get_fdata()

def angular_error(PEa, PEb, halfPi=True):
    # PEa = preprocessing.normalize(PEa, norm='l1', axis=1)
    # PEb = preprocessing.normalize(PEb, norm='l1', axis=1)
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
            np.square(PEa[..., 1] - PEb[..., 1]) + \
            np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)

    ang = 2 * np.real(np.arcsin(chord/2))

    if halfPi:
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]
    return np.degrees(ang)

err = angular_error(a,b)
ape =  np.abs(err)

#mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_mask.nii.gz'))
mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/tractsegout/TOM_trackings/*_mask.nii.gz'))
values = []
for mask_file in mask_files:
    mask = nib.load(mask_file).get_fdata()
    roi = err * mask
    mean_roi = np.nanmean(roi)
    #print(mask_file, mean_roi)
    values.append(mean_roi)

print(values)
print(min(values))
print(max(values))
