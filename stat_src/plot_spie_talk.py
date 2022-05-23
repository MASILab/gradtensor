import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap()

true_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()
corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_study/Lest_md.nii').get_fdata()
ape_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/src/md_ape_image.nii').get_fdata()

plt.figure();
slice_idx = 45
slice = true_md[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,1)
plt.axis('off')
#cmap = plt.get_cmap('parula')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=parula)
plt.title('Intended FA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

slice = corpt_md[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=parula)
plt.title('Achieved FA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

slice = ape_md[slice_idx,:,:]
m = 0
M = 0.05
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,3,3)
plt.axis('off')
cmap = plt.get_cmap('bwr')
cmap = cmap.reversed() 
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Non-linear effects')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
plt.show()
