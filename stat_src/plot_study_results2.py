import numpy as np
import nibabel as nib
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap

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



plt.figure();
slice_idx = 45;
slice = true_fa[:,slice_idx,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("Ground Truth FA", fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corpt_change = abs(corpt_fa - true_fa);
slice = corpt_change[:,slice_idx,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 FA \n Corrupt - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Corruption and Ground Truth")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

sm_corr_change = abs(sm_corr_fa-true_fa);
slice = sm_corr_change[:,slice_idx,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,4)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

bx_corr_change = abs(bx_corr_fa-true_fa);
slice = bx_corr_change[:,slice_idx,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,3)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 FA \n Empirical - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Empirical Bammer Correction \n Method and Ground Truth")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)





slice = true_md[:,slice_idx,:]
m = 0
M = 0.003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,5)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("Ground Truth MD", fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.set_yticklabels([0.0, 0.5, 1, 1.5, 2, 2.5, 3])
a.ax.tick_params(labelsize=12)
a.set_label(u"\u03bcm2/s", size = 12)

corpt_change = abs(corpt_md - true_md);
slice = corpt_change[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,6)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 MD \n Corrupt - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs MD change in \n Corruption and Ground Truth", fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.set_yticklabels([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
a.ax.tick_params(labelsize=12)
a.set_label(u"\u03bcm2/s", size = 12)

sm_corr_change = abs(sm_corr_md-true_md);
slice = sm_corr_change[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,8)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 MD \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs MD change in \n Shortcut Correction Method \n and Ground Truth", fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.set_yticklabels([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
a.ax.tick_params(labelsize=12)
a.set_label(u"\u03bcm2/s", size = 12)

bx_corr_change = abs(bx_corr_md-true_md);
slice = bx_corr_change[:,slice_idx,:]
m = 0
M = 0.0003
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 MD \n Empirical - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs MD change in \n Empirical Bammer Correction \n Method and Ground Truth")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.set_yticklabels([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
a.ax.tick_params(labelsize=12)
a.set_label(u"\u03bcm2/s", size = 12)



def angular_error_real(PEa, PEb, halfPi=True):
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

slice = true_pev[:, slice_idx,:]
m = 0
M = 1
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,9)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title('Ground Truth PEV', fontdict = {'fontsize' : 15})
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('degrees', size = 12)

corpt_change = angular_error_real(corpt_pev,true_pev);
slice = corpt_change[:, slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,10)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 PEV \n Corrupt - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Corruption and Ground Truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('degrees', size = 12)

sm_corr_change = angular_error_real(sm_corr_pev,true_pev);
slice = sm_corr_change[:, slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,12)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 PEV \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Shortcut Simple Correction Method \n and Ground Truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('degrees', size = 12)

bx_corr_change = angular_error_real(bx_corr_pev,true_pev);
slice = bx_corr_change[:, slice_idx,:]
m = 0
M = 5
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap= parula)
plt.title("\u0394 PEV \n Empirical - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title('Angular change in \n Empirical Bammer Correction \n Method  and Ground Truth')
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)
a.ax.tick_params(labelsize=12)
a.set_label('degrees', size = 12)

plt.subplots_adjust(left=0.01, bottom=0.05, right=0.92, top=0.9, wspace=0.3, hspace=0.2)
plt.show()
