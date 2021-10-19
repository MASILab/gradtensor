import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt


MD_sm = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_est/sm_md.nii').get_fdata()
MD_inv_sm = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_est/2_inv_sm__md.nii').get_fdata()
