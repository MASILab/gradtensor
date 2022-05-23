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
corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_fa.nii').get_fdata()
noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_Nest_fa.nii').get_fdata()
lcorpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_Lest_fa.nii').get_fdata()
mask = nib.load('/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii').get_fdata()

inf_corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/uncorrected_fa.nii').get_fdata()
inf_noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/uncorrected_Nest_fa.nii').get_fdata()

sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_corrected_fa.nii').get_fdata()
Nsm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_Ncorrected_fa.nii').get_fdata()
bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/emp/emp_corrected_fa.nii').get_fdata()
Nbx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/Nemp/Nemp_corrected_fa.nii').get_fdata()
#Nbx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/test_noise/test_fa.nii').get_fdata()
#Nbx_corr_fa_inf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/Nemp/Nemp_corrected_fa.nii').get_fdata()

def plot_one_snr(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa =  (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    df_ape_fa = df_ape_fa.assign(x = x)
    return df_ape_fa

inf_noise = plot_one_snr(inf_noise_fa,true_fa,'Uncorrected SNR=inf','inf_noise')
inf_lr_noise = plot_one_snr(inf_corpt_fa,true_fa,'Uncorrected SNR=inf','inf_lr_noise')
inf_lr =  plot_one_snr(inf_corpt_fa,inf_noise_fa,'Uncorrected SNR=inf','inf_lr')

noise = plot_one_snr(noise_fa,true_fa,'Uncorrected SNR=100','noise')
lr_noise = plot_one_snr(corpt_fa,true_fa,'Uncorrected SNR=100','lr_noise')
lr =  plot_one_snr(corpt_fa,noise_fa,'Uncorrected SNR=100','lr')

#corr_noise = plot_one_snr(Nbx_corr_fa,true_fa,'Approx corrected SNR=30','corr_noise')
#corr_lr_noise = plot_one_snr(bx_corr_fa,true_fa, 'Approx corrected SNR=30','corr_lr_noise')
#corr_lr = plot_one_snr(bx_corr_fa,noise_fa,'Approx corrected SNR=30','corr_lr')

corr_noise = plot_one_snr(Nsm_corr_fa,true_fa,'Approx corrected SNR=100','corr_noise')
corr_lr_noise = plot_one_snr(sm_corr_fa,true_fa, 'Approx corrected SNR=100','corr_lr_noise')
corr_lr = plot_one_snr(sm_corr_fa,noise_fa,'Approx corrected SNR=100','corr_lr')

plt.figure()
sns.set(font_scale = 1.3)
pd_data = pd.concat([inf_noise, noise, corr_noise, inf_lr_noise, lr_noise,corr_lr_noise,inf_lr, lr, corr_lr])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
ax = sns.boxplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge =False,width=.5)
ax.set(xlabel=' ', ylabel = 'Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )')
ax.set_xticklabels(['','noise FA - GT FA','','','noise+LR FA - GT FA','','','noise+LR FA - noise FA',''])
plt.legend(loc='upper right')
plt.title('Effect of noise and LR on FA in whole brain')
plt.ylim([-5,150])
plt.show()
