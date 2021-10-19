import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

MD_Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_sm/Lest_corr_sm_md.nii').get_fdata()
MD_true = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii').get_fdata()

plt.figure()
slice_idx = 45

slice = MD_Lest[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of Signal after Simple Correction')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

slice = MD_true[slice_idx,:,:]
m = slice.min()
M = slice.max()
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD of ground truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

diff_uncorr = MD_Lest - MD_true
slice = diff_uncorr[slice_idx,:,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD Simple Corrected - Ground truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

plt.show()

