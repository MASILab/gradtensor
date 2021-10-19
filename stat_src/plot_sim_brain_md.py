import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

# MD
#Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/posB_Lest_md.nii').get_fdata()
#uncorr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/brain_posB_md.nii').get_fdata()

#est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/posB_est_md.nii').get_fdata()
#corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii').get_fdata()

#corr_md_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/posB_corr_md.nii').get_fdata()

est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/_OUTPUTS_signal_estimate/est_md.nii').get_fdata()
uncorr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/brainA__fa.nii').get_fdata()

Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/_OUTPUTS_signal_estimate/co_Lest_og_ob_md.nii').get_fdata()
corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()

#m = 0
#M = 0.0033
m = 0
M = 0.003
slice_idx = 45 
slice = uncorr[:,slice_idx,:]
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
a = plt.figure()
plt.subplot(3,3,1)
plt.axis('off') 
cmap = plt.get_cmap('viridis') 
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = est[:, slice_idx,:]
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = corr[:, slice_idx,:] 
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3)) 
slice = np.nan_to_num(slice)
plt.subplot(3,3,4)
plt.axis('off')
cmap = plt.get_cmap('viridis') 
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)


slice = Lest[:, slice_idx,:] 
m = slice.min()
print(m)
M = slice.max()
print(M)
slice = np.flip(np.rot90(slice,3)) 
slice = np.nan_to_num(slice)
plt.subplot(3,3,5)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

#m = -0.04
#M = 0.1
diff1 = uncorr - est
slice = diff1[:, slice_idx,:]
slice = np.flip(np.rot90(slice,3)) 
slice = np.nan_to_num(slice)
plt.subplot(3,3,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

diff2 = corr - Lest
slice = diff2[:, slice_idx,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3)) 
slice = np.nan_to_num(slice)
plt.subplot(3,3,6)
plt.axis('off')
cmap = plt.get_cmap('viridis') 
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

#m = 0
#M = 0.0004
diff3= uncorr - corr
slice = diff3[:, slice_idx,:]
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,7)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

diff4 =  Lest - est
m = 0
M = 0.0001
slice = diff4[:, slice_idx,:]
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,8)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

a.suptitle('FA')
#a.text(0.1, 0.5, 'w Lr                                    w/o Lr', ha='center', va='center', rotation='vertical')
a.text(0.45, 0.09, '      DWI                        Sim DWI                    Difference', ha='center', va='center')
a.subplots_adjust(right=0.8,hspace=0.02,wspace=0)
cbar_ax = a.add_axes([0.85, 0.15, 0.02, 0.7])
a.colorbar(im, cbar_ax)
plt.show()
