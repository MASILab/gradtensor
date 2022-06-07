import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.io
import h5py
import mat73

SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNRinf.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR30.mat')

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

def snr_brain(FA_sim_noise_x,FA_true_x,label,x):
    err_fa = FA_sim_noise_x - FA_true_x
    pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(pe_fa)
    cube = np.reshape(ape_fa, [50, 50, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[2,]
    df_brain = pd.DataFrame(brain).assign(x = x)
    df_brain = df_brain.assign(label = label)
    median = df_brain.groupby(['label', 'x'])[0].median().values
    return df_brain, median

def pev_brain(FA_sim_noise_x,FA_true_x,label,x):
    err_fa = angular_error(FA_sim_noise_x,FA_true_x)
    #pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(err_fa)
    cube = np.reshape(ape_fa, [50, 50, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[2,]
    df_brain = pd.DataFrame(brain).assign(x = x)
    df_brain = df_brain.assign(label = label)
    median = df_brain.groupby(['label', 'x'])[0].median().values
    return df_brain, median

def violin_plot(metric):
    if metric != 'PEV':
        noise_braininf,minf_noise= snr_brain(SNRinf[metric+'_sim_noise_x'],SNRinf[metric+'_true_x'],'Uncorrected SNR=inf','inf_noise');
        noiseLR_braininf,minf_lr_noise = snr_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected SNR=inf','inf_lr_noise');
        LR_braininf,minf_lr = snr_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_sim_noise_x'],'Uncorrected SNR=inf','inf_lr');

        noise_brain30,mnoise = snr_brain(SNR30[metric+'_sim_noise_x'],SNR30[metric+'_true_x'],'Uncorrected SNR=30','noise');
        noiseLR_brain30,mlr_noise = snr_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected SNR=30','lr_noise');
        LR_brain30,mlr = snr_brain(SNRinf[metric+'_sim_corpt_x'],SNR30[metric+'_sim_noise_x'],'Uncorrected SNR=30','lr');

        corr_noise_brain30,mcorr_noise = snr_brain(SNR30[metric+'_corr_noise_bx_x'],SNR30[metric+'_true_x'],'Emp corrected SNR=30','corr_noise');
        corr_noiseLR_brain30,mcorr_lr_noise = snr_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_true_x'],'Emp corrected SNR=30','corr_lr_noise');
        corr_LR_brain30,mcorr_lr = snr_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp corrected SNR=30','corr_lr');
        pd_data = pd.concat([noise_braininf,noise_brain30,corr_noise_brain30,noiseLR_braininf,noiseLR_brain30,corr_noiseLR_brain30,LR_braininf,LR_brain30,corr_LR_brain30])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
        m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
        mL1 = [(np.round(s, 2)) for s in m1]

    else:
        noise_braininf,minf_noise= pev_brain(SNRinf[metric+'_sim_noise_x'],SNRinf[metric+'_true_x'],'Uncorrected SNR=inf','inf_noise');
        noiseLR_braininf,minf_lr_noise = pev_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected SNR=inf','inf_lr_noise');
        LR_braininf,minf_lr = pev_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_sim_noise_x'],'Uncorrected SNR=inf','inf_lr');

        noise_brain30,mnoise = pev_brain(SNR30[metric+'_sim_noise_x'],SNR30[metric+'_true_x'],'Uncorrected SNR=30','noise');
        noiseLR_brain30,mlr_noise = pev_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected SNR=30','lr_noise');
        LR_brain30,mlr = pev_brain(SNRinf[metric+'_sim_corpt_x'],SNR30[metric+'_sim_noise_x'],'Uncorrected SNR=30','lr');

        corr_noise_brain30,mcorr_noise = pev_brain(SNR30[metric+'_corr_noise_bx_x'],SNR30[metric+'_true_x'],'Emp corrected SNR=30','corr_noise');
        corr_noiseLR_brain30,mcorr_lr_noise = pev_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_true_x'],'Emp corrected SNR=30','corr_lr_noise');
        corr_LR_brain30,mcorr_lr = pev_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp corrected SNR=30','corr_lr');
        pd_data = pd.concat([noise_braininf,noise_brain30,corr_noise_brain30,noiseLR_braininf,noiseLR_brain30,corr_noiseLR_brain30,LR_braininf,LR_brain30,corr_LR_brain30])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
        m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
        mL1 = [(np.round(s, 2)) for s in m1]
    return pd_data, mL1

plt.subplots(3,1,figsize=(15,20))
sns.set(font_scale = 1.3)
sns.set_style("white")
palette = {'Uncorrected SNR=inf': 'crimson', 'Uncorrected SNR=30': 'cornflowerblue', 'Emp corrected SNR=30': 'limegreen'}
plt.subplot(3,1,1)
fa_pd_data,mL1 = violin_plot('FA')
ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=1500,palette=palette)
ax.set_ylim([-10,100])
ax.legend(loc='upper right')
ax.set(xlabel=' ', ylabel = 'APE in FA (%)')
for xtick in ax.get_xticks():
    ax.text(xtick-.2,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise FA - GT FA','','','noise+LR FA - GT FA','','','noise+LR FA - noise FA',''])
plt.grid()

plt.subplot(3,1,2)
md_pd_data,mL1 = violin_plot('MD')
ax = sns.violinplot(data=md_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=10000,palette=palette)
ax.set_ylim([-15,100])
ax.get_legend().remove()
ax.set(xlabel=' ', ylabel = 'APE in MD (%)')
for xtick in ax.get_xticks():
    ax.text(xtick-.2,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise MD - GT MD','','','noise+LR MD - GT MD','','','noise+LR MD - noise MD',''])
plt.grid()

plt.subplot(3,1,3)
v1_pd_data,mL1 = violin_plot('PEV')
ax = sns.violinplot(data=v1_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=1500,palette=palette)
ax.set_ylim([-10,100])
ax.set(xlabel=' ', ylabel = 'APE in V1 (degrees)')
ax.get_legend().remove()

for xtick in ax.get_xticks():
    ax.text(xtick-.2,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
            horizontalalignment='center',size='x-small',color='k',weight='semibold')
ax.set_xticklabels(['','noise PEV - GT PEV','','','noise+LR PEV - GT PEV','','','noise+LR PEV - noise PEV',''])
plt.grid()  

plt.tight_layout()
plt.show()
