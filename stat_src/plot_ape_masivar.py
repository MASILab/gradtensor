import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


mask = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/mask.nii').get_fdata()

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

def pev_df(corpt_fa,true_fa,label,x):
    err_fa = angular_error(corpt_fa,true_fa)
    #pe_fa =  (err_fa / true_fa) * 100
    ape_fa = (err_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    df_ape_fa = df_ape_fa.assign(x = x)
    median = df_ape_fa.groupby(['label', 'x'])[0].median().values
    return df_ape_fa, median


def violin_plot(metric):
    true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+metric+'.nii').get_fdata()
    #true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__'+metric+'.nii').get_fdata()
    corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    lcorpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Lest_'+metric+'.nii').get_fdata()

    inf_corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    inf_noise_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()

    sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_corrected_'+metric+'.nii').get_fdata()
    Nsm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_Ncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    Nbx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/Nemp/Nemp_corrected_'+metric+'.nii').get_fdata()
    if metric != 'primary_eigvec':
        inf_noise,minf_noise = plot_one_snr(inf_noise_fa,true_fa,'Uncorrected SNR=inf','inf_noise')
        inf_lr_noise,minf_lr_noise = plot_one_snr(inf_corpt_fa,true_fa,'Uncorrected SNR=inf','inf_lr_noise')
        inf_lr,minf_lr =  plot_one_snr(inf_corpt_fa,inf_noise_fa,'Uncorrected SNR=inf','inf_lr')

        noise,mnoise = plot_one_snr(noise_fa,true_fa,'Uncorrected SNR=30','noise')
        lr_noise,mlr_noise = plot_one_snr(corpt_fa,true_fa,'Uncorrected SNR=30','lr_noise')
        lr,mlr =  plot_one_snr(corpt_fa,noise_fa,'Uncorrected SNR=30','lr')

        corr_noise,mcorr_noise= plot_one_snr(Nbx_corr_fa,true_fa,'Emp corrected SNR=30','corr_noise')
        corr_lr_noise,mcorr_lr_noise = plot_one_snr(bx_corr_fa,true_fa, 'Emp corrected SNR=30','corr_lr_noise')
        corr_lr,mcorr_lr = plot_one_snr(bx_corr_fa,noise_fa,'Emp corrected SNR=30','corr_lr')
        pd_data = pd.concat([inf_noise, noise, corr_noise, inf_lr_noise, lr_noise,corr_lr_noise,inf_lr, lr, corr_lr])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
        m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
        mL1 = [(np.round(s, 2)) for s in m1]

    else:
        inf_noise,minf_noise = pev_df(inf_noise_fa,true_fa,'Uncorrected SNR=inf','inf_noise')
        inf_lr_noise,minf_lr_noise = pev_df(inf_corpt_fa,true_fa,'Uncorrected SNR=inf','inf_lr_noise')
        inf_lr,minf_lr = pev_df(inf_corpt_fa,inf_noise_fa,'Uncorrected SNR=inf','inf_lr')

        noise,mnoise = pev_df(noise_fa,true_fa,'Uncorrected SNR=30','noise')
        lr_noise,mlr_noise = pev_df(corpt_fa,true_fa,'Uncorrected SNR=30','lr_noise')
        lr,mlr =  pev_df(corpt_fa,noise_fa,'Uncorrected SNR=30','lr')

        corr_noise,mcorr_noise= pev_df(Nbx_corr_fa,true_fa,'Emp corrected SNR=30','corr_noise')
        corr_lr_noise,mcorr_lr_noise = pev_df(bx_corr_fa,true_fa, 'Emp corrected SNR=30','corr_lr_noise')
        corr_lr,mcorr_lr = pev_df(bx_corr_fa,noise_fa,'Emp corrected SNR=30','corr_lr')
        pd_data = pd.concat([inf_noise, noise, corr_noise, inf_lr_noise, lr_noise,corr_lr_noise,inf_lr, lr, corr_lr])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
        m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
        mL1 = [(np.round(s, 2)) for s in m1]
    return pd_data, mL1

plt.subplots(3,1,figsize=(15,20))
sns.set(font_scale = 2)
sns.set_style("white")
palette = {'Uncorrected SNR=inf': 'crimson', 'Uncorrected SNR=30': 'cornflowerblue', 'Emp corrected SNR=30': 'limegreen'}
plt.subplot(3,1,1)
fa_pd_data,mL1 = violin_plot('fa')
ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=2000,palette=palette,inner="box",bw=0.08)
ax.set_ylim([-10,100])
ax.legend(loc='upper right', prop={'size': 17})
ax.set(xlabel=' ', ylabel = 'APE in FA (%)')
#for xtick in ax.get_xticks():
#    ax.text(xtick-.3,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
#            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise FA - GT FA','','','noise+LR FA - GT FA','','','noise+LR FA - noise FA',''])
plt.grid()

plt.subplot(3,1,2)
md_pd_data,mL1 = violin_plot('md')
ax = sns.violinplot(data=md_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=10000,palette=palette,inner="box",bw=0.08)
ax.set_ylim([-15,100])
ax.get_legend().remove()
ax.set(xlabel=' ', ylabel = 'APE in MD (%)')
#for xtick in ax.get_xticks():
#    ax.text(xtick-.4,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
#            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise MD - GT MD','','','noise+LR MD - GT MD','','','noise+LR MD - noise MD',''])
plt.grid()

plt.subplot(3,1,3)
v1_pd_data,mL1 = violin_plot('primary_eigvec')
ax = sns.violinplot(data=v1_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=1500,palette=palette,inner="box",bw=0.08)
ax.set_ylim([-10,100])
ax.set(xlabel=' ', ylabel = 'AE in V1 (°)')
ax.get_legend().remove()

#for xtick in ax.get_xticks():
#    ax.text(xtick-.3,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
#            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise V1 - GT V1','','','noise+LR V1 - GT V1','','','noise+LR V1 - noise V1',''])

plt.grid()  
plt.tight_layout()
plt.show()
