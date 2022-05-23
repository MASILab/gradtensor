import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap

true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__fa.nii').get_fdata()
corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_1/uncorrected_fa.nii').get_fdata()
sm_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/real_corr_sm_fa.nii').get_fdata()
bx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/real_corr_bx_fa.nii').get_fdata()

true_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__md.nii').get_fdata()
corpt_md = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_1/uncorrected_md.nii').get_fdata()
sm_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/real_corr_sm_md.nii').get_fdata()
bx_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/real_corr_bx_md.nii').get_fdata()

true_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__primary_eigvec.nii').get_fdata()
corpt_pev = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_1/uncorrected_primary_eigvec.nii').get_fdata()
sm_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/real_corr_sm_primary_eigvec.nii').get_fdata()
bx_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/real_corr_bx_primary_eigvec.nii').get_fdata()

mask =  nib.load('/home/local/VANDERBILT/kanakap/mask.nii').get_fdata()
img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__fa.nii')

true_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__ad.nii').get_fdata()
corpt_ad = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_1/uncorrected_ad.nii').get_fdata()
sm_corr_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/real_corr_sm_ad.nii').get_fdata()
bx_corr_ad = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/real_corr_bx_ad.nii').get_fdata()

true_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__rd.nii').get_fdata()
corpt_rd = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_1/uncorrected_rd.nii').get_fdata()
sm_corr_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/real_corr_sm_rd.nii').get_fdata()
bx_corr_rd = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/real_corr_bx_rd.nii').get_fdata()

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

fig,faxs = plt.subplots(3,5);
plt.subplots_adjust(wspace= 0, hspace= 0.1, bottom=0.5, right=0.95, left=0.1)
slice_idx = 52
mask = mask[slice_idx,:,:]
mask = np.array(mask, dtype=bool)

err_fa = corpt_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = true_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
im = faxs[0,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
faxs[0,0].axis('off')
cb = fig.colorbar(im, ax = faxs[0, 0],shrink=0.9,aspect=5,ticks=[0.0,1])
cb.ax.tick_params(labelsize=15)

err_md = corpt_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = true_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,1].axis('off')
im = faxs[0,1].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

err_ad = corpt_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = true_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,2].axis('off')
faxs[0,2].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

err_rd = corpt_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = true_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,3].axis('off')
im = faxs[0,3].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[0, 1:4],shrink=0.9,aspect=5,ticks=[0.0, 0.003])
cb.ax.tick_params(labelsize=15)

pev_err = angular_error(corpt_pev, true_pev)
slice = pev_err[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[0,4].axis('off')
im = faxs[0,4].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[0, 4],shrink=0.9,aspect=5,ticks=[0.0,5.0])
cb.ax.tick_params(labelsize=15)

err_fa = bx_corr_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = bx_corr_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,0].axis('off')
im = faxs[1,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[1, 0],shrink=0.9,aspect=5,ticks=[0.0,1])
cb.ax.tick_params(labelsize=15)

err_md = bx_corr_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = bx_corr_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M =  0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,1].axis('off')
faxs[1,1].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

pev_err = angular_error(bx_corr_pev, true_pev)
slice = pev_err[slice_idx,:,:]
#mk = np.array([mask,mask,mask])
#idx = [1, 2, 0]
#idx = np.empty_like(permutation)
#idx[permutation] = np.arange(len(permutation))
#mk1 = mk[:, idx]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,4].axis('off')
im = faxs[1,4].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[1, 4],shrink=0.9,aspect=5,ticks=[0.0,5.0])
cb.ax.tick_params(labelsize=15)

err_ad = bx_corr_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = bx_corr_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,2].axis('off')
faxs[1,2].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

err_rd = bx_corr_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = bx_corr_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[1,3].axis('off')
im = faxs[1,3].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[1, 1:4],shrink=0.9,aspect=5,ticks=[0.0, 0.003])
cb.ax.tick_params(labelsize=15)

err_fa = sm_corr_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
slice = sm_corr_fa[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,0].axis('off')
im = faxs[2,0].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[2, 0],shrink=0.9,aspect=5,ticks=[0.0,1])
cb.ax.tick_params(labelsize=15)

err_md = sm_corr_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)
slice = sm_corr_md[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,1].axis('off')
faxs[2,1].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

err_ad = sm_corr_ad - true_ad
pe_ad = 100 * err_ad / true_ad
ape_ad = np.abs(pe_ad)
slice = sm_corr_ad[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,2].axis('off')
faxs[2,2].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')

err_rd = sm_corr_rd - true_rd
pe_rd = 100 * err_rd / true_rd
ape_rd = np.abs(pe_rd)
slice = sm_corr_rd[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,3].axis('off')
im = faxs[2,3].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[2, 1:4],shrink=0.9,aspect=5,ticks=[0.0, 0.003])
cb.ax.tick_params(labelsize=15)

pev_err = angular_error(sm_corr_pev, true_pev)
slice = pev_err[slice_idx,:,:]
slice = np.where(mask, slice, 0)
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
faxs[2,4].axis('off')
im = faxs[2,4].imshow(slice, vmin=m, vmax=M, cmap= 'viridis')
cb = fig.colorbar(im, ax = faxs[2, 4],shrink=0.9,aspect=5,ticks=[0.0,5.0])
cb.ax.tick_params(labelsize=15)
plt.show()



