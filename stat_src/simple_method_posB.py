import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# FA
uncorr_md_posA = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posA_md.nii').get_fdata()
uncorr_md_posB = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posB_md.nii').get_fdata()

corr_md_posA = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_future_fieldmap/posA_corr_md.nii').get_fdata()
corr_md_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/posB_corr_md.nii').get_fdata()


uncorr_fa_posA = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posA_fa.nii').get_fdata()
uncorr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/INPUTS/brain_posB_fa.nii').get_fdata()

corr_fa_posA = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_future_fieldmap/posA_corr_fa.nii').get_fdata()
corr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/posB_corr_fa.nii').get_fdata()
Lcorr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_est/sm_fa.nii').get_fdata()

corr_sig = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/new_meth_corrected_sig.nii').get_fdata()
uncorr_sig = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.nii').get_fdata()

#for analysis
rot_uncorr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/rot_posB_fa.nii').get_fdata()
rot_corr_fa_posB =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/rot_posB_fa.nii').get_fdata()

rot45_uncorr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/rot45_posB_fa.nii').get_fdata()
rot45_corr_fa_posB = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/rot45_posB_fa.nii').get_fdata()

#for i in range(96):
#    for j in range(96):
#        for k in range(68):
#            if corr_fa_posB[i,j,k] > 0.8:
#                print(corr_fa_posB[i,j,k])
#                print(i,j,k)
#print(uncorr_fa_posB[42,38,35])
#print(Lcorr_fa_posB[42,38,35])
print(rot_uncorr_fa_posB[42,38,35])
print(rot_corr_fa_posB[42,38,35])
print(']]]')
print(rot45_uncorr_fa_posB[42,38,35])
print(rot45_corr_fa_posB[42,38,35])

exit()
slice_idx = 45
plt.figure()

uncorr = uncorr_sig[slice_idx,:,:15]
slice = uncorr[slice_idx,:,:]
m = 0
M = 2.5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(4,4,1)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Uncorrected DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

corr = corr_sig[:,:,:,15]
slice = corr[slice_idx,:,:]
m = 0
M = 2.5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(4,4,9)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Corrected DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_sig = corr_sig - uncorr_sig
diff_sig = diff_sig[:,:,:,15]
slice = diff_sig[slice_idx,:,:]
m = 0
M = 0.35
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(4,4,13)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Corrected - Uncorrected sig')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = uncorr_md_posB[slice_idx,:,:] 
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3)) 
plt.subplot(4,4,3)
plt.axis('off')
cmap = plt.get_cmap('viridis') 
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD uncorrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


slice = corr_md_posB[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
plt.subplot(4,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD corrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_md = corr_md_posB - uncorr_md_posB
m = slice.min()
M = slice.max()
slice = diff_md[slice_idx,:,:]
slice = np.flip(np.rot90(slice,3))
plt.subplot(4,4,12)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD corrected - uncorrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


slice = uncorr_fa_posB[slice_idx,:,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
plt.subplot(4,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA uncorrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


slice = corr_fa_posB[slice_idx,:,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
plt.subplot(4,4,15)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA corrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


diff_fa = corr_fa_posB - uncorr_fa_posB
m = -0.04
M = 0.04
slice = diff_fa[slice_idx,:,:]
slice = np.flip(np.rot90(slice,3))
plt.subplot(4,4,16)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('FA corrected - uncorrected posB')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

plt.show()
