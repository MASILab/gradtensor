import nibabel as nib
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#for analysis
true =nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_fa.nii').get_fdata()

rot_corr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/rot_Lest_corr_fa.nii').get_fdata()
rot_uncorr_fa_posB =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot_Lest_fa.nii').get_fdata()

rot45_corr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/rot45_Lest_corr_fa.nii').get_fdata()
rot45_uncorr_fa_posB =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot45_Lest_fa.nii').get_fdata()
rot_corr_bx_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_rot_bax_method/rot_Lest_corr_bax_fa.nii').get_fdata()

rot45_corr_bx_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_rot_bax_method/rot45_Lest_corr_bax_fa.nii').get_fdata()


#print(rot45_uncorr_fa_posB[73,53,43])
#print(rot45_corr_bx_fa_posB[73,53,43])
#exit()
print('true')
print(true[73,53,43])
#print('simple')
#for i in range(96):
#    for j in range(96):
#        for k in range(68):
#            #print(math.isnan(rot_uncorr_fa_posB[i,j,k]))
#            #print(math.isnan(rot_corr_fa_posB[i,j,k]))
#            if not (math.isnan(rot_uncorr_fa_posB[i,j,k]) or math.isnan(rot_corr_bx_fa_posB[i,j,k])):
#                if round(rot_uncorr_fa_posB[i,j,k],1) != round(rot_corr_bx_fa_posB[i,j,k],1):
#                    print(i,j,k)
#                    print(round(rot_uncorr_fa_posB[i,j,k],1))
#                    print(round(rot_corr_bx_fa_posB[i,j,k],1))
#exit()
# 42,38,35
print(rot_uncorr_fa_posB[73,53,43])
print(rot_corr_fa_posB[73,53,43])
print('simple45')
print(rot45_uncorr_fa_posB[73,53,43])
print(rot45_corr_fa_posB[73,53,43])


rot_corr_bx_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_rot_bax_method/rot_Lest_corr_bax_fa.nii').get_fdata()

rot45_corr_bx_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_rot_bax_method/rot45_Lest_corr_bax_fa.nii').get_fdata()

print('bx')
print(rot_uncorr_fa_posB[73,53,43])
print(rot_corr_bx_fa_posB[73,53,43])
print('bx45')
print(rot45_uncorr_fa_posB[73,53,43])
print(rot45_corr_bx_fa_posB[73,53,43])
