import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


MD_sm = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_est/sm_md.nii').get_fdata()
MD_sm_inv = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_simple_method_est/2_sm_md.nii').get_fdata()

plt.figure()
slice_idx = 45

diff_uncorr = MD_sm - MD_sm_inv 
slice = diff_uncorr[slice_idx,:,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(1,1,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('MD(Est D^ uncorr) - MD(D)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=ax)
plt.show()
