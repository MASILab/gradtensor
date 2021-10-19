import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

ten_raw = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/3_dwmri_md.nii').get_fdata()
est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_posA_3tb_est_sig.nii').get_fdata()
ten_est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_est_md.nii').get_fdata()

ten_L = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()
Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_posA_3tb_Lest_sig.nii').get_fdata()
ten_Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_Lest_md.nii').get_fdata()

plt.figure()
slice_idx = 45
slice = ten_raw[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of raw DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

est = est[:,:,:,15]
slice = est[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,2)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Simulate DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


slice = ten_est[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of Simulate DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_est = ten_raw - ten_est
slice = diff_est[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,4)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD Difference')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_L[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,5)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of DWI with Lr')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

Lest = Lest[:,:,:,15]
slice = Lest[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,6)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Simulate DWI with Lr')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_Lest[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of Simulate DWI with Lr')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_Lest = ten_L - ten_Lest
slice = diff_Lest[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,8)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD difference')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff = ten_L - ten_raw
slice = diff[slice_idx,:,:]
m = 0
M = 0.0001
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,9)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD difference')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

#diff_sig = est - Lest
diff_sig = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_posA_3tb_diff_sig.nii').get_fdata()
diff_sig = diff_sig[:,:,:,15]
slice = diff_sig[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,10)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Sim sig difference')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_L = ten_Lest - ten_est
slice = diff_L[slice_idx,:,:]
m = 0
M = 0.0001
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD difference')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

plt.show()

