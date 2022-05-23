import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import operator as op
import xml.etree.cElementTree as et
import pandas as pd
import scipy.io as sio

corpt = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/FAdiff_avg_wmlabels_ape.nii.gz').get_fdata()

sm = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/FAdiff_avg_wmlabels_ape_SMCORR.nii.gz').get_fdata()

bx = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/FAdiff_avg_wmlabels_ape_BXCORR.nii.gz').get_fdata()


fig,faxs = plt.subplots(1,3);
slice_idx = 52
plt.subplots_adjust(wspace= 0.1, hspace= 0.1, bottom=0.5, right=0.95, left=0.1)
slice = corpt[slice_idx,:,:]
m = 0
M = 2
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
faxs[0].axis('off')

slice = bx[slice_idx,:,:]
m = 0
M = 2
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1].axis('off')
faxs[1].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

slice = sm[slice_idx,:,:]
m = 0
M = 2
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2].axis('off')
im = faxs[2].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[ 0:3],shrink=0.5,aspect=20,ticks=[0.0,2.0])
cb.ax.tick_params(labelsize=13)

plt.show()
