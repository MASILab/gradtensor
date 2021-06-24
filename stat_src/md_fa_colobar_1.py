import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

# FA
uncorr_md_posA = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_FA.nii.gz').get_fdata()
uncorr_md_posB = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/INPUTS/dti_FA.nii.gz').get_fdata()
uncorr_md_posC = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/INPUTS/dti_FA.nii.gz').get_fdata()

corr_md_posA = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()
corr_md_posB = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_fa.nii').get_fdata()
corr_md_posC = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/OUTPUTS_future_fieldmap/p_3tb_posC_mask_fa.nii').get_fdata()


slice_idx = 50 #maybe pick a middle slice
slice = uncorr_md_posA[:,slice_idx,:]# get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
a = plt.figure()
plt.subplot(3,3,1)
plt.axis('off')
m = -0.01
M = 0.1 #this min and max needs to be set for each metric
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = uncorr_md_posB[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(3,3,4)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = uncorr_md_posC[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(3,3,7)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = corr_md_posA[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)


slice = corr_md_posB[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,5)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

slice = corr_md_posC[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,8)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

m = -0.02
M = 0.2
diff_posA = uncorr_md_posA - corr_md_posA
slice = diff_posA[:, slice_idx,:]
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,3)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

diff_posB = uncorr_md_posB - corr_md_posB
slice = diff_posB[:, slice_idx,:]
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,6)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)

diff_posC = uncorr_md_posC - corr_md_posC
slice = diff_posC[:, slice_idx,:]
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
slice = np.nan_to_num(slice)
plt.subplot(3,3,9)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)


#cbar_ax = a.add_axes([0.85, 0.15, 0.05, 0.7])
a.suptitle('FA')
a.text(0.1, 0.5, 'posC                 posB                     posA', ha='center', va='center', rotation='vertical')
a.text(0.45, 0.07, 'Uncorrected              Corrected            Difference', ha='center', va='center')
a.subplots_adjust(right=0.8,hspace=0.02,wspace=0)
cbar_ax = a.add_axes([0.85, 0.15, 0.02, 0.7])
a.colorbar(im, cbar_ax)
plt.show()
