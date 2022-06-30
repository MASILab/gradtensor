import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import nmmn.plots
#parula=nmmn.plots.parulacmap() # for MATLAB's cmap
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__md.nii').get_fdata()
corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_md.nii').get_fdata()
noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Nest_md.nii').get_fdata()
lcorpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Lest_md.nii').get_fdata()
#mask = nib.load('/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii').get_fdata()
mask = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/mask.nii').get_fdata()

inf_corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_md.nii').get_fdata()
inf_noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_Nest_md.nii').get_fdata()

sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_corrected_md.nii').get_fdata()
Nsm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_Ncorrected_md.nii').get_fdata()
bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/emp/emp_corrected_md.nii').get_fdata()
Nbx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/Nemp/Nemp_corrected_md.nii').get_fdata()
#Nbx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/test_noise/test_fa.nii').get_fdata()
#Nbx_corr_fa_inf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1/Nemp/Nemp_corrected_fa.nii').get_fdata()


def plot_one_snr(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa =  (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [i for i in ape_fa if math.isnan(i) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(x = x)
    df_ape_fa = df_ape_fa.assign(label = label)
    median = df_ape_fa.groupby(['label', 'x'])[0].median().values
    return df_ape_fa, median

inf_noise,minf_noise = plot_one_snr(inf_noise_fa,true_fa,'Uncorrected SNR=inf','inf_noise')
inf_lr_noise,minf_lr_noise = plot_one_snr(inf_corpt_fa,true_fa,'Uncorrected SNR=inf','inf_lr_noise')
inf_lr,minf_lr =  plot_one_snr(inf_corpt_fa,inf_noise_fa,'Uncorrected SNR=inf','inf_lr')

noise,mnoise = plot_one_snr(noise_fa,true_fa,'Uncorrected SNR=30','noise')
lr_noise,mlr_noise = plot_one_snr(corpt_fa,true_fa,'Uncorrected SNR=30','lr_noise')
lr,mlr =  plot_one_snr(corpt_fa,noise_fa,'Uncorrected SNR=30','lr')

corr_noise,mcorr_noise= plot_one_snr(Nbx_corr_fa,true_fa,'Emp corrected SNR=30','corr_noise')
corr_lr_noise,mcorr_lr_noise = plot_one_snr(bx_corr_fa,true_fa, 'Emp corrected SNR=30','corr_lr_noise')
corr_lr,mcorr_lr = plot_one_snr(bx_corr_fa,noise_fa,'Emp corrected SNR=30','corr_lr')
plt.figure(1)
sns.set(font_scale = 1.3)
pd_data = pd.concat([inf_noise, noise, corr_noise, inf_lr_noise, lr_noise,corr_lr_noise,inf_lr, lr, corr_lr])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
sns.set_style("white")
palette = {'Uncorrected SNR=inf': 'crimson', 'Uncorrected SNR=30': 'cornflowerblue', 'Emp corrected SNR=30': 'limegreen'}
plt.subplot(3,1,2)
ax = sns.violinplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=10000,palette=palette)
ax.set(xlabel=' ', ylabel = 'APE in MD (%)')
ax.set_xticklabels(['','noise FA - GT FA','','','noise+LR FA - GT FA','','','noise+LR FA - noise FA',''])
m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
print(m1)
mL1 = [(np.round(s, 2)) for s in m1]
print(mL1)
for xtick in ax.get_xticks():
    ax.text(xtick-.2,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

#plt.legend(loc='upper right')
#plt.title('Effect of noise and LR on FA in whole brain')
ax.get_legend().remove()
plt.ylim([-15,100])
plt.grid()
plt.tight_layout()
plt.show()
