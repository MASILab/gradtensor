import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()
corpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_study/Lest_fa.nii').get_fdata()
sm_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_study/Lest_corr_sm_fa.nii').get_fdata()
bx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_study/Lest_corr_bx_fa.nii').get_fdata()

true_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()
corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_study/Lest_md.nii').get_fdata()
sm_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_study/Lest_corr_sm_md.nii').get_fdata()
bx_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_study/Lest_corr_bx_md.nii').get_fdata()

true_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_primary_eigvec.nii').get_fdata()
corpt_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_study/Lest_primary_eigvec.nii').get_fdata()
sm_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_study/Lest_corr_sm_primary_eigvec.nii').get_fdata()
bx_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_study/Lest_corr_bx_primary_eigvec.nii').get_fdata()
img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii')

err_fa = corpt_fa - true_fa
valid = ~np.isnan(err_fa) * ~np.isnan(true_fa) * true_fa!=0
true_fa = true_fa[valid]
err_fa = err_fa[valid]
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)

ape_img_fa = nib.Nifti1Image(ape_fa, img.affine, img.header)
nib.save(ape_img_fa, 'fa_ape_image.nii')

err_md = corpt_md - true_md
pe_md = 100 * err_md / true_md
ape_md = np.abs(pe_md)

ape_img_md = nib.Nifti1Image(ape_md, img.affine, img.header)
nib.save(ape_img_md, 'md_ape_image.nii')


"""
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

plt.figure();
slice_idx = 45
slice = true_fa[:,slice_idx,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Ground Truth FA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

diff_corpt_fa = abs(corpt_fa - true_fa)
slice = diff_corpt_fa[:,slice_idx,:]
m =  slice.min()
M =  0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Corrupt FA and Ground Truth FA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)

diff_bx_fa = abs(bx_corr_fa - true_fa)
slice = diff_bx_fa[:,slice_idx,:]
m =  slice.min()
M =  0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Empirical Bammer \n Correction Method FA and Ground truth FA ')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)

diff_sm_fa = abs(sm_corr_fa - true_fa)
slice = diff_sm_fa[:,slice_idx,:]
m =  slice.min()
M =  0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,4)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Shortcut Method \n Correction Method FA and Ground truth FA')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)


slice = true_md[slice_idx,:,:]
m = slice.min()
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,5)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Ground Truth MD')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

diff_corpt_md = abs(corpt_md-true_md)
slice = diff_corpt_md[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,6)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Corrpution MD and Ground truth MD')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

diff_corpt_md = abs(corpt_md-true_md)
slice = diff_corpt_md[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Empirical Bammer Correction \n Method MD and Ground truth MD')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)

diff_sm_md = abs(sm_corr_md-true_md)
slice = diff_sm_md[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,8)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Abs Change in Shortcut Method  Correction \n MD and Ground truth MD')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('mm2/s', size = 12)


slice_idx = 45;
slice = true_pev[:,slice_idx,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,9)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Grouth Truth PEV')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

diff_corpt_pev = angular_error(true_pev,corpt_pev)
slice = diff_corpt_pev[:,slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,10)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Angular Change in Corrput PEV \n and Ground Truth PEV')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

diff_bx_pev = angular_error(true_pev,bx_corr_pev,halfPi=True)
slice = diff_bx_pev[:,slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Angular Change in Empirical Bammer \n Correction Method and Ground truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

diff_sm_pev = angular_error(true_pev,sm_corr_pev,halfPi=True)
slice = diff_sm_pev[:,slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,12)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
plt.title('Angular Change in Shortcut Correction \n Method and Ground truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

plt.show()
"""
