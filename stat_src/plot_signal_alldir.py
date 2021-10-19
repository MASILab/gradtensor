import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

corr_sig = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap/new_meth_corrected_sig.nii').get_fdata()


fig, axs = plt.subplots(5,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
slice_idx = 45

for i in range(25):

    corr = corr_sig[:,:,:,i]
    slice = corr[slice_idx,:,:]
    m = 0
    M = 2.5
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    axs[i].axis('off')
    cmap = plt.get_cmap('gray')
    axs[i].imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
    #axs[i].set_title(str(250+i))

plt.show()
