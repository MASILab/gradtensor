import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap

corpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_fa.nii').get_fdata()
emp_fa =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_emp_fa.nii').get_fdata()
appr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_appr_fa.nii').get_fdata()

corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_md.nii').get_fdata()
emp_md =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_emp_md.nii').get_fdata()
appr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_appr_md.nii').get_fdata()

corpt_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_rd.nii').get_fdata()
emp_rd =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_emp_rd.nii').get_fdata()
appr_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_appr_rd.nii').get_fdata()

corpt_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_ad.nii').get_fdata()
emp_ad =  nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_emp_ad.nii').get_fdata()
appr_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor/stat_src/cohend_image_snr30_appr_ad.nii').get_fdata()


fig,faxs = plt.subplots(3,4);
plt.subplots_adjust(wspace= 0, hspace= 0.1, bottom=0.5, right=0.95, left=0.1)
slice_idx = 52

slice = corpt_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,0].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
faxs[0,0].axis('off')

slice = corpt_md[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,1].axis('off')
faxs[0,1].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = corpt_ad[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,2].axis('off')
faxs[0,2].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = corpt_rd[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,3].axis('off')
im = faxs[0,3].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
cb = fig.colorbar(im, ax = faxs[0, :4],shrink=0.9,aspect=5,ticks=[-3,0.0,3])
cb.ax.tick_params(labelsize=15)

slice = emp_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,0].axis('off')
faxs[1,0].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = emp_md[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,1].axis('off')
faxs[1,1].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = emp_ad[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,2].axis('off')
faxs[1,2].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = emp_rd[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,3].axis('off')
im = faxs[1,3].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
cb = fig.colorbar(im, ax = faxs[1, :4],shrink=0.9,aspect=5,ticks=[-3,0.0,3])
cb.ax.tick_params(labelsize=15)

slice = appr_fa[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,0].axis('off')
faxs[2,0].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = appr_md[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,1].axis('off')
faxs[2,1].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = appr_ad[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,2].axis('off')
faxs[2,2].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')

slice = appr_rd[slice_idx,:,:]
m = -3
M = 3
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,3].axis('off')
im = faxs[2,3].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
cb = fig.colorbar(im, ax = faxs[2, :4],shrink=0.9,aspect=5,ticks=[-3,0.0,3])
cb.ax.tick_params(labelsize=15)

plt.show()



