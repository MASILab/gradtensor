import nibabel as nib
import numpy as np
import scipy.io as sio


true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()

corpt_fa_dict = {}
for i in range(0,4):
    pre = '/nfs/masi/kanakap/projects/LR/aggregate_study/'
    post = '/uncorrected_fa.nii'
    outdir = pre + 'OUTPUT_' +sr(i) + post
    corpt_fa = nib.load(outdir)
    corpt_fa_dict[i] = corpt_fa


