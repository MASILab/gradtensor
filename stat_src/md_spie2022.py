import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

ten_raw = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/brain_posB_md.nii').get_fdata()
est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/posB_3tb_est_sig.nii').get_fdata()
ten_est = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/est_md.nii').get_fdata()

ten_L = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii').get_fdata()
Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/posB_3tb_Lest_sig.nii').get_fdata()
ten_Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/Lest_md.nii').get_fdata()
MD_est_d_uncorr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/Lest_og_ob_md.nii').get_fdata()
MD_sm_d_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_est/sm_md.nii').get_fdata()

MD_fy_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_anisotropy_method/Lest_corr_md.nii').get_fdata()

MD_bl_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_bl_method/Lest_corr_md.nii').get_fdata()

fig = plt.figure()
slice_idx = 45

slice = MD_est_d_uncorr[slice_idx,:,:]
m = slice.min()
M = 0.003
print(M)
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
sb = plt.subplot(1,3,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
#plt.title('MD (Est D^uncorr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=23)
a.set_label('mm2/s', size = 25)

slice = ten_L[slice_idx,:,:]
m = slice.min()
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
sb = plt.subplot(1,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
#plt.title('MD (D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=23)
a.set_label('mm2/s', size = 25)

diff_uncorr = MD_est_d_uncorr - ten_L
slice = diff_uncorr[slice_idx,:,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
sb = plt.subplot(1,3,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
img = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(img, cax=ax)
a.ax.tick_params(labelsize=23)
a.set_label('mm2/s', size = 25)
#img = plt.gca()
#plt.title('MD(Est D^ uncorr) - MD(D)')
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#fig.subplots_adjust(right=0.8,hspace=0.02,wspace=0)
#cbar_ax = fig.add_axes([0.75, 0.15, 0.02, 0.2])
#a = fig.colorbar(im, cbar_ax)
#a.set_label('mm2/s')
#cbar = plt.colorbar(im)
#tick_font_size = 12
#cbar_ax.tick_params(labelsize=tick_font_size)
#fig.tight_layout()
plt.show()
