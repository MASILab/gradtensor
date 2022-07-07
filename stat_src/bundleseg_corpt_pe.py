import os
import sys
from glob import glob
import nibabel as nib
import numpy as np


a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_est/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

def angular_error(PEa, PEb, halfPi=True):
    # PEa = preprocessing.normalize(PEa, norm='l1', axis=1)
    # PEb = preprocessing.normalize(PEb, norm='l1', axis=1)
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
            np.square(PEa[..., 1] - PEb[..., 1]) + \
            np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)

    ang = 2 * np.real(np.arcsin(chord/2))

    if halfPi:
        print(ang)
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]
        print(ang)
    return np.degrees(ang)

err = angular_error(a,b)
ape =  np.abs(err)

mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_mask.nii.gz'))

for mask_file in mask_files:
    mask = nib.load(mask_file).get_fdata()
    roi = ape * mask
    mean_roi = np.nanmean(roi)
    print(mask_file, mean_roi)
