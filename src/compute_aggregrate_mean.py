import nibabel as nib
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap
from sklearn.metrics import mean_squared_error

def rmse_voxelwise(imga, imgb):
    rmse_img = np.zeros((112,112,54))
    for x in range(112):
        for y in range(112):
            for z in range(54):
                err = np.mean((imga[x,y,z] - imgb[x,y,z]) ** 2)
                err = np.sqrt(err)
                rmse_img[x,y,z] = err
    return rmse_img

def mse_voxelwise(imga, imgb):
    mse_img = np.zeros((112,112,54))
    for x in range(112):
        for y in range(112):
            for z in range(54):
                err = np.mean((imga[x,y,z] - imgb[x,y,z]) ** 2)
                err = np.sqrt(err)
                mse_img[x,y,z] = err
    return mse_img

mse_img_dict = {}
def varience_voxelwise(mean_img, img_dict):
    for i in range(1,11):
        mse_img = mse_voxelwise(mean_img, img_dict[i])
        mse_img_dict[i] = mse_img
    mse_img_list = np.array(list(mse_img_dict.values()))
    varience_voxelwise = np.nanmean(mse_img_list,0)
    return varience_voxelwise

mask =  nib.load('/home/local/VANDERBILT/kanakap/mask.nii').get_fdata()
#true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()
true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__fa.nii').get_fdata()
true_fa = np.where(mask, true_fa, 0)

corpt_fa_dict = {}
cpre = '/nfs/masi/kanakap/projects/LR/aggregate_study/'
cpost = '/uncorrected_fa.nii'
for i in range(1,11):
    outdir = cpre + 'OUTPUT_masivar_SNR100_' +str(i) + cpost
    corpt_fa = nib.load(outdir).get_fdata()
    corpt_fa = np.where(mask, corpt_fa, 0)
    corpt_fa_dict[i] = corpt_fa

corpt_fa_list = np.array(list(corpt_fa_dict.values()))
corpt_fa_mean = np.nanmean(corpt_fa_list,0)
corpt_fa_variance = varience_voxelwise(corpt_fa_mean, corpt_fa_dict) #np.nanvar(corpt_fa_list,0) #/ corpt_fa_mean

corr_sm_fa_dict = {}
pre = '/nfs/masi/kanakap/projects/LR/aggregate_study/'
post = '/approx_corrected_fa.nii'
for i in range(1,11):
    outdir = pre + 'OUTPUT_masivar_SNR100_' +str(i) + post
    corr_fa = nib.load(outdir).get_fdata()
    corr_fa = np.where(mask, corr_fa, 0)
    corr_sm_fa_dict[i] = corr_fa

corr_sm_fa_list = np.array(list(corr_sm_fa_dict.values()))
corr_sm_fa_mean = np.nanmean(corr_sm_fa_list,0)
corr_sm_fa_variance = varience_voxelwise(corr_sm_fa_mean, corr_sm_fa_dict)
#corr_sm_fa_variance = np.nanvar(corr_sm_fa_list,0) #/ corr_sm_fa_mean

epre = '/nfs/masi/kanakap/projects/LR/aggregate_study/'
epost = '/emp_corrected_fa.nii'
corr_bx_fa_dict = {}
for i in range(1,11):
    outdir = epre + 'OUTPUT_masivar_SNR100_' +str(i) + epost
    corr_fa = nib.load(outdir).get_fdata()
    corr_fa = np.where(mask, corr_fa, 0)
    corr_bx_fa_dict[i] = corr_fa

corr_bx_fa_list = np.array(list(corr_bx_fa_dict.values()))
corr_bx_fa_mean = np.nanmean(corr_bx_fa_list,0)
corr_bx_fa_variance = varience_voxelwise(corr_bx_fa_mean,corr_bx_fa_dict)
#corr_bx_fa_variance = np.nanvar(corr_bx_fa_list,0) #/ corr_bx_fa_mean
mask =  nib.load('/home/local/VANDERBILT/kanakap/mask.nii').get_fdata()
plt.figure();
slice_idx = 52;

corpt_bias = rmse_voxelwise(corpt_fa_mean,true_fa)
slice = corpt_bias[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,1)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("Ground Truth FA", fontdict = {'fontsize' : 15})
plt.title("Corrpution Bias")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

slice = corpt_fa_variance[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,2)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Corrupt - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Corruption and Ground Truth")
plt.title("Corruption Varience")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corr_bx_bias = rmse_voxelwise(corr_bx_fa_mean,true_fa)
slice = corr_bx_bias[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,5)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap='viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("Bias Emp. Correction")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

slice = corr_bx_fa_variance[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,6)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("Varience Emp. Correction")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corr_bx_fa_diff = corr_bx_bias - corpt_bias
slice = corr_bx_fa_diff[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,7)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("\u0394 Bais \n Emp. Correction \n and Corrupt")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corr_bx_fa_diff = corr_bx_fa_variance - corpt_fa_variance
slice = corr_bx_fa_diff[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,8)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("\u0394 Varience \n Emp. Correction \n and Corrupt")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corr_sm_bais = rmse_voxelwise(corr_sm_fa_mean,true_fa)
slice = corr_sm_bais[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,9)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("Approx. Correction Bais")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

slice = corr_sm_fa_variance[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,10)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
#plt.title("Abs FA change in \n Shortcut Correction Method \n and Ground Truth")
plt.title("Approx. Correction Varience")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)


corr_sm_fa_diff = corr_sm_bais - corpt_bias
slice = corr_sm_fa_diff[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,11)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
plt.title("\u0394 Bais \n Approx. Correction \n and Corrupt")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

corr_sm_fa_diff = (corr_sm_fa_variance - corpt_fa_variance)
slice = corr_sm_fa_diff[slice_idx,:,:]
m = 0
M = 0.04
slice = np.flip(np.rot90(slice,3))
slice = np.nan_to_num(slice)
plt.subplot(3,4,12)
plt.axis('off')
cmap = plt.get_cmap('viridis')
im = plt.imshow((slice), vmin=m, vmax=M, cmap= 'viridis')
#plt.title("\u0394 FA \n Approx. - Ground Truth", fontdict = {'fontsize' : 15})
plt.title("\u0394 Varience \n Approx. Correction\n and Corrupt")
divider = make_axes_locatable(plt.gca())
ax = divider.append_axes("right", size="5%", pad=0.05)
a = plt.colorbar(im, cax=ax)

plt.subplots_adjust(left=0.01, bottom=0.05, right=0.92, top=0.9, wspace=0.3, hspace=0.2)
plt.show()
