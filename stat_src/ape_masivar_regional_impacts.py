import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


metrics = ['fa','md','v1','ad','rd']
for i in metrics:
     locals()['{0}'.format(i)] = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/'+ i + '_diff_gm_seg.nii.gz').get_fdata()

metrics = [fa,md,ad,rd,v1]
fig,faxs = plt.subplots(2,5);
#plt.subplots_adjust(wspace= 0, hspace= 0.1, bottom=0.5, right=0.95, left=0.1)
slice_idx = 40
for i in range(len(metrics)):
    a  = metrics[i]
    slice = a[slice_idx,:,:]
    m = 0
    M = 10
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    faxs[0,i].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
    faxs[0,i].axis('off')

metrics = ['fa','md','v1','ad','rd']
for i in metrics:
     locals()['{0}'.format(i)] = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/'+ i + '_diff_wm_seg.nii.gz').get_fdata()

metrics = [fa,md,ad,rd,v1]
slice_idx = 40
for i in range(len(metrics)):
    a  = metrics[i]
    slice = a[slice_idx,:,:]
    m = 0
    M = 3
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    faxs[1,i].imshow(slice, vmin=m, vmax=M, cmap= 'bwr')
    faxs[1,i].axis('off')
plt.show()
