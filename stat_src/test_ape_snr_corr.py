import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


noise_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_Nest_fa.nii').get_fdata()
bx_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_corrected_fa.nii').get_fdata()

noise_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/uncorrected_Nest_fa.nii').get_fdata()
bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/approx_corrected_fa.nii').get_fdata()

noise_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/uncorrected_Nest_fa.nii').get_fdata()
bx_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/approx_corrected_fa.nii').get_fdata()
mask = nib.load('/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii').get_fdata()


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


corr_lr100 = plot_one_snr(bx_corr_fa100,noise_fa100,'SNR=100','emp_corr_lr_100')
corr_lr30 = plot_one_snr(bx_corr_fa30,noise_fa30,'SNR=30','emp_corr_lr_30')
corr_lr10 = plot_one_snr(bx_corr_fa10,noise_fa10,'SNR=10','emp_corr_lr_10')


acorr_lr100 = plot_one_snr(sm_corr_fa100,noise_fa100,'SNR=100','appr_corr_lr_100')
acorr_lr30 = plot_one_snr(sm_corr_fa30,noise_fa30,'SNR=30','appr_corr_lr_30')
acorr_lr10 = plot_one_snr(sm_corr_fa10,noise_fa10,'SNR=10','appr_corr_lr_10')

plt.figure()
sns.set(font_scale = 1.3)
pd_data = pd.concat([corr_lr100, corr_lr30, corr_lr10, acorr_lr100, acorr_lr30, acorr_lr10])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
ax = sns.boxplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge =False,width=.5)
ax.set(xlabel=' ', ylabel = 'Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )')
ax.set_xticklabels(['','Emp. Correction','','','Approximate Correction',''])
plt.legend(loc='upper right')
plt.title('Effect of noise and LR on FA in whole brain')
plt.ylim([-0.5,2])
plt.show()
