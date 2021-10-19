import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

MD_est_d_corr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri_md.nii').get_fdata()
MD_est_d_uncorr = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/Lest_aug15_md.nii').get_fdata()

plt.figure()
slice_idx = 45
diff = MD_est_d_uncorr - MD_est_d_corr
slice = diff[slice_idx,:,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
#plt.subplot(1,3,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
#plt.title('FA (Est D^uncorr)')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
cb = plt.colorbar(im, cax=ax)
cb.set_label('mm2/s')
plt.show()
#plt.savefig('fafig1.png')

