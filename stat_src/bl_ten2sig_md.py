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
MD_est_d_uncorr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_Lest_og_ob_md.nii').get_fdata()
MD_sm_d_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_est/2_inv_sm_md.nii').get_fdata()

MD_fy_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_anisotropy_method/Lest_corr_md.nii').get_fdata()


plt.figure()
slice_idx = 45

#simulated w no L
est = est[:,:,:,15]
slice = est[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,1)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Obs Ideal DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_est[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_raw[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff = ten_est - ten_L
slice = diff[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,4)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est^D) - MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

#simulated w L
Lest = Lest[:,:,:,15]
slice = Lest[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,5)
plt.axis('off')
cmap = plt.get_cmap('gray')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Obs Physical DWI')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)


slice = MD_est_d_uncorr[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,6)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^uncorr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_raw[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_uncorr = MD_est_d_uncorr - ten_L
slice = diff_uncorr[slice_idx,:,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,8)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD(Est D^ uncorr) - MD(D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_Lest[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,10)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^corr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_L[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_corr = ten_Lest - ten_L
slice = diff_corr[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,12)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD(Est D^corr) - MD(D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = MD_sm_d_corr[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,14)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^sm_corr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_L[slice_idx,:,:]
m = slice.max()
M = slice.min()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,15)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_sm_corr = MD_sm_d_corr - ten_L
slice = diff_sm_corr[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,16)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^sm_corr) - MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = MD_fy_corr[slice_idx,:,:]
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
m = 0
M = 0.003
plt.subplot(5,4,18)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^fy_corr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

slice = ten_L[slice_idx,:,:]
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,19)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

diff_fy_corr = MD_fy_corr - ten_L
slice = diff_fy_corr[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(5,4,20)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD (Est D^fy_corr) - MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)

#diff in signal
diff_sig = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/3_posA_3tb_diff_sig.nii').get_fdata()
diff_sig = diff_sig[:,:,:,15]
slice = diff_sig[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(4,4,9)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Obs Ideal - Obs Phy sig')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)
plt.show()

