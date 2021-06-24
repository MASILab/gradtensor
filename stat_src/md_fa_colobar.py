import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

# MD
uncorr_md_posA = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_MD.nii.gz').get_fdata()
uncorr_md_posB = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/INPUTS/dti_MD.nii.gz').get_fdata()
uncorr_md_posC = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/INPUTS/dti_MD.nii.gz').get_fdata()

corr_md_posA = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()
corr_md_posB = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii').get_fdata()
corr_md_posC = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/OUTPUTS_future_fieldmap/p_3tb_posC_mask_md.nii').get_fdata()

fig, axes = plt.subplots(2,3, figsize=(8,8))

slice_idx = 32 #maybe pick a middle slice
slice = uncorr_md_posA[:, :, slice_idx]# get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
m = 0
M = 0.0001 #this min and max needs to be set for each metric
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
img1 = axes[0,0].imshow(slice, vmin=m, vmax=M, cmap=cmap)


slice = uncorr_md_posB[:, :, slice_idx] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
img2 = axes[0,1].imshow(slice, vmin=m, vmax=M, cmap=cmap)

slice = uncorr_md_posC[:, :, slice_idx] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
img3 = axes[1,0].imshow(slice, vmin=m, vmax=M, cmap=cmap)

slice = corr_md_posA[:, :, slice_idx] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
img4 = axes[1,1].imshow(slice, vmin=m, vmax=M, cmap=cmap)


slice = corr_md_posB[:, :, slice_idx] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
img5 = axes[0,0].imshow(slice, vmin=m, vmax=M, cmap=cmap)

slice = corr_md_posC[:, :, slice_idx] # get axial slice
slice = np.flip(np.rot90(slice)) #nibabel load the nifti weird so gotta flip and rotate
img6 = axes[1,1].imshow(slice, vmin=m, vmax=M, cmap=cmap)

fig.colorbar(img1, ax=axes[:, 1])
plt.show()
