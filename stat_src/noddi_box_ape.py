import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import sys

uncorr_output_dir = sys.argv[1]
corr_output_dir = sys.argv[2]

uncorr_ficvf = nib.load(uncorr_output_dir + '/noddi_ficvf.nii').get_fdata()
uncorr_fiso = nib.load(uncorr_output_dir + '/noddi_fiso.nii').get_fdata()
uncorr_odi = nib.load(uncorr_output_dir + '/noddi_odi.nii').get_fdata()

corr_ficvf = nib.load(corr_output_dir + '/noddi_ficvf.nii').get_fdata()
corr_fiso = nib.load(corr_output_dir + '/noddi_fiso.nii').get_fdata()
corr_odi = nib.load(corr_output_dir + '/noddi_odi.nii').get_fdata()

def calc_ape(corpt_fa,true_fa):
    mask = np.logical_and(corpt_fa > 0, true_fa > 0)
    err_fa = corpt_fa - true_fa
    pe_fa =  (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = ape_fa[mask]
    #ape_fa = ape_fa.reshape(-1)
    return ape_fa

def calc_diff(corpt_fa,true_fa):
    mask = np.logical_and(corpt_fa > 0, true_fa > 0)
    err_fa = np.abs(corpt_fa - true_fa)
    err_fa[np.isnan(err_fa)] = 0
    err_fa = err_fa[mask]
    #ape_fa = ape_fa.reshape(-1)
    return err_fa

def calc_vrmse(a,b):
    mask = np.logical_and(a > 0, b > 0)
    RMSE_all = np.zeros((a.shape[0],a.shape[1],a.shape[2]))
    for x in range(a.shape[0]):
        for y in range(a.shape[1]):
            for z in range(a.shape[2]):
                MSE = np.square(b[x,y,z] - a[x,y,z]).mean()
                RMSE = math.sqrt(MSE)
                RMSE_all[x,y,z] = RMSE
    RMSE_all = RMSE_all[mask]
    return RMSE_all

def voxelwise_ape(a,b):
    ape = []
    for x in range(a.shape[0]):
        for y in range(a.shape[1]):
            for z in range(a.shape[2]):
                err_fa = (b[x,y,z] - a[x,y,z])
                pe_fa = ( err_fa / a[x,y,z]) * 100
                ape_fa = np.abs(pe_fa)
                ape.append(ape_fa)
    return ape

def plot_ape(image,i):
    slice = image[slice_idx, :,:]
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    plt.subplot(3,3,i)
    plt.axis('off')
    cmap = plt.get_cmap('viridis')
    m = 0
    M = 100
    im = plt.imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
    divider = make_axes_locatable(plt.gca())
    ax = divider.append_axes("right", size="5%", pad=0.05)
    a = plt.colorbar(im, cax=ax)

ape_ficvf = calc_ape(uncorr_ficvf,corr_ficvf)
ape_fisco = calc_ape(uncorr_fiso,corr_fiso)
ape_odi = calc_ape(uncorr_odi,corr_odi)
data = [ape_ficvf, ape_fisco, ape_odi]

plt.rcParams.update({'font.size':20})
plt.boxplot(data,flierprops={'marker': '.', 'markersize': 2})
plt.ylim([-2, 100])
plt.xticks([1, 2, 3], ['iVF', 'cVF', 'ODI'])
plt.show()
