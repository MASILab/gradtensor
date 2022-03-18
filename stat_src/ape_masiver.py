import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap

true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__fa.nii').get_fdata()
corpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/masiver_fa.nii').get_fdata()
sm_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/approx/Lest_corr_sm_fa.nii').get_fdata()
bx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/full/Lest_corr_bx_fa.nii').get_fdata()

true_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__md.nii').get_fdata()
corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/masiver_md.nii').get_fdata()
sm_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/approx/Lest_corr_sm_md.nii').get_fdata()
bx_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/full/Lest_corr_bx_md.nii').get_fdata()

true_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__primary_eigvec.nii').get_fdata()
corpt_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/masiver_primary_eigvec.nii').get_fdata()
sm_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/approx/Lest_corr_sm_primary_eigvec.nii').get_fdata()
bx_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/full/Lest_corr_bx_primary_eigvec.nii').get_fdata()
mask =  nib.load('/home/local/VANDERBILT/kanakap/mask.nii').get_fdata()
img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__fa.nii')

true_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__ad.nii').get_fdata()
corpt_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/masiver_ad.nii').get_fdata()
sm_corr_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/approx/Lest_corr_sm_ad.nii').get_fdata()
bx_corr_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/full/Lest_corr_bx_ad.nii').get_fdata()

true_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__rd.nii').get_fdata()
corpt_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/masiver_rd.nii').get_fdata()
sm_corr_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/approx/Lest_corr_sm_rd.nii').get_fdata()
bx_corr_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/full/Lest_corr_bx_rd.nii').get_fdata()

fig,fax=plt.subplots(3,5);
slice_idx = 52
mask = mask[slice_idx,:,:]
mask = np.array(mask, dtype=bool)

err_fa = corpt_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = ape_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,1)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution FA", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_md = corpt_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = ape_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,2)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution MD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

def angular_error(PEa, PEb, halfPi=True):
    # PEa = preprocessing.normalize(PEa, norm='l1', axis=1)
    # PEb = preprocessing.normalize(PEb, norm='l1', axis=1)
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
            np.square(PEa[..., 1] - PEb[..., 1]) + \
            np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)

    ang = 2 * np.real(np.arcsin(chord/2))

    if halfPi:
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]

    return np.degrees(ang)

pev_err = angular_error(corpt_pev, true_pev)
slice = pev_err[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,5)
plt.axis('off')
#cmap = plt.get_cmap('hot')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution PEV", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Corruption and Ground Truth')
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)
#a.ax.tick_params(labelsize=12)
#a.set_label('degrees', size = 12)

err_ad = corpt_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = ape_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,3)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution AD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_rd = corpt_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = ape_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,4)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution RD", fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = fig.colorbar(im, cax=ax,shrink=0.6)

err_fa = corpt_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = ape_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,1)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Corrpution FA", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)


err_fa = bx_corr_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = ape_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,6)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Emp. Correction FA", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_md = bx_corr_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = ape_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,7)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Emp. Correction MD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

pev_err = angular_error(bx_corr_pev, true_pev)
slice = pev_err[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,10)
plt.axis('off')
#cmap = plt.get_cmap('hot')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Emp. Correction PEV", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Corruption and Ground Truth')
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)
#a.ax.tick_params(labelsize=12)
#a.set_label('degrees', size = 12)

err_ad = bx_corr_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = ape_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,8)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Emp. Correction AD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_rd = sm_corr_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = ape_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,9)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction RD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_fa = sm_corr_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = ape_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,11)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction FA", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

err_md = sm_corr_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = ape_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,12)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction MD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)

pev_err = angular_error(sm_corr_pev, true_pev)
slice = pev_err[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,15)
plt.axis('off')
#cmap = plt.get_cmap('hot')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction PEV", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Corruption and Ground Truth')
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)
#a.ax.tick_params(labelsize=12)
#a.set_label('degrees', size = 12)

err_ad = sm_corr_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = ape_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,13)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction AD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(, cax=ax)

err_rd = sm_corr_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = ape_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 10
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,5,14)
plt.axis('off')
im = plt.imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Approx. Correction RD", fontdict = {'fontsize' : 15})
#divider = make_axes_locatable(plt.gca())
#ax = divider.append_axes("right", size="5%", pad=0.05)
#a = plt.colorbar(im, cax=ax)
plt.subplots_adjust(hspace=0, wspace=0,bottom=0.39,left=0.026,right=0.96)
plt.show()



