import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# FA
uncorr_md_posA = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posA_fa.nii').get_fdata()
uncorr_md_posB = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posB_fa.nii').get_fdata()

corr_md_posA = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_future_fieldmap/posA_corr_fa.nii').get_fdata()
corr_md_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/posB_corr_fa.nii').get_fdata()

slice_idx = 45 #maybe pick a middle slice
slice = uncorr_md_posA[:,slice_idx,:]# get axial slice
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
a = plt.figure()
plt.subplot(2,3,1)
plt.axis('off')
m = 0
M = 1 #this min and max needs to be set for each metric
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA uncorrected posA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = uncorr_md_posB[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(2,3,4)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA uncorrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = corr_md_posA[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(2,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA corrected posA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = corr_md_posB[:, slice_idx,:] # get axial slice
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(2,3,5)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA corrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_posA = uncorr_md_posA - corr_md_posA
slice = diff_posA[:, slice_idx,:]
m = -0.04
M = 0.04
print(m)
print(M)
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(2,3,3)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA uncorrected - FA corrected posA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


diff_posB = uncorr_md_posB - corr_md_posB
slice = diff_posB[:, slice_idx,:]
m = -0.04
M = 0.04
print(m)
print(M)
slice = np.flip(np.rot90(slice,3)) #nibabel load the nifti weird so gotta flip and rotate
plt.subplot(2,3,6)
plt.axis('off')
cmap = plt.get_cmap('viridis') #this is default but you can choose whatever
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA uncorrected - FA corrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

a.suptitle('FA')
#a.text(0.1, 0.5, 'posB                                    posA', ha='center', va='center', rotation='vertical')
#a.text(0.45, 0.09, 'Uncorrected              Corrected            Difference', ha='center', va='center')
#a.subplots_adjust(right=0.8,hspace=0.02,wspace=0)
#cbar_ax = a.add_axes([0.85, 0.15, 0.02, 0.7])
#a.colorbar(im, cbar_ax)
plt.show()
