import os
import sys
import subprocess
from glob import glob
import multiprocessing as mp

import nibabel as nib
import numpy as np

from dipy.tracking.utils import density_map


#trk_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_cleaned.trk'))
#trk_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/tractsegout/TOM_trackings/*.trk'))
trk_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs052/ses-s1Bx1/tracto_op_lr_corr_1/bundles/*_cleaned.trk'))

for trk_file in trk_files:

        print('Converting {}'.format(trk_file))
        trk = nib.streamlines.load(trk_file)

        tdi = density_map(streamlines=trk.streamlines, affine=trk.affine, vol_dims=trk.header['dimensions'])
        tdi_nii = nib.Nifti1Image(tdi, trk.affine)
        tdi_file = trk_file.replace('.trk', '_tdi.nii.gz')
        nib.save(tdi_nii, tdi_file)

        thres = 0.05*np.percentile(tdi[tdi > 0], 99)
        mask = (tdi > thres).astype(int)
        mask_file = tdi_file.replace('.nii.gz', '_mask.nii.gz')
        mask_nii = nib.Nifti1Image(mask, trk.affine)
        nib.save(mask_nii, mask_file)
