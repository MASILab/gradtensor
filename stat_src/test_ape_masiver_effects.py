import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap

corpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr10.nii').get_fdata()
emp_fa =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr100_emp_fa.nii').get_fdata()
appr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr10_appr_fa.nii').get_fdata()

corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr10_md.nii').get_fdata()

corpt_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr10_rd.nii').get_fdata()

fig,faxs = plt.subplots(3,4);
plt.subplots_adjust(wspace= 0, hspace= 0.1, bottom=0.5, right=0.95, left=0.1)
slice_idx = 52

slice = corpt_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
faxs[0,0].axis('off')

slice = corpt_md[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,1].axis('off')
faxs[0,1].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')


slice = corpt_rd[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,3].axis('off')
im = faxs[0,3].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[0, :4],shrink=0.9,aspect=5,ticks=[0.0,10.0])
cb.ax.tick_params(labelsize=15)

slice = emp_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,0].axis('off')
faxs[1,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')


slice = appr_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,0].axis('off')
faxs[2,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')


plt.show()



