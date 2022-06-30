import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import nmmn.plots
#parula=nmmn.plots.parulacmap() # for MATLAB's cmap
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math



true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__fa.nii').get_fdata()
corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_fa.nii').get_fdata()
noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
lcorpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_Lest_fa.nii').get_fdata()
mask = nib.load('/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii').get_fdata()

inf_corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/uncorrected_fa.nii').get_fdata()
inf_noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/uncorrected_Nest_fa.nii').get_fdata()

sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/approx_corrected_fa.nii').get_fdata()
Nsm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/approx_Ncorrected_fa.nii').get_fdata()
bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
Nbx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/Nemp/Nemp_corrected_fa.nii').get_fdata()

def plot_one_snr(corpt_fa,true_fa,label,x):
    #corpt_fa[np.isnan(corpt_fa)] = 0
    #true_fa[np.isnan(true_fa)] = 0
    err_fa = corpt_fa - true_fa
    pe_fa =  (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    m = np.median(ape_fa)
    return ape_fa,m


noise,mnoise = plot_one_snr(noise_fa,true_fa,'Uncorrected SNR=30','noise')
lr_noise,mlr_noise = plot_one_snr(corpt_fa,true_fa,'Uncorrected SNR=30','lr_noise')
lr,mlr =  plot_one_snr(corpt_fa,noise_fa,'Uncorrected SNR=30','lr')

tips = sns.load_dataset("tips")

def plot_comparison(x, title):
    fig, ax = plt.subplots(3, 1, sharex=True)
    sns.distplot(x, ax=ax[0])
    ax[0].set_title('Histogram + KDE')
    sns.boxplot(x, ax=ax[1])
    ax[1].set_title('Boxplot')
    sns.violinplot(x, ax=ax[2],gridsize=1000)
    ax[2].set_title('Violin plot')
    fig.suptitle(title, fontsize=16)
    ax[2].set_xlabel('Abs percent error in FA')
    plt.rcParams.update({'font.size': 40})
    plt.show()

plot_comparison(noise,'Effect of noise on FA at SNR = 30')
#plot_comparison(tips,'lala')
#N = 10
#sample_gaussian = np.random.normal(size=N)
#plot_comparison(sample_gaussian, 'Standard Normal Distribution')

"""
print(mnoise)
plt.figure()
sns.set_theme(style="whitegrid")
sns.violinplot(np.array(lr_noise))
#plt.xlim([0, 150])
plt.show()
"""
