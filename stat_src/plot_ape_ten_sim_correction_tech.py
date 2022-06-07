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
SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR30.mat')
SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR10.mat')

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

def split_voilin(delta=0.02,inner="quart"):
    for ii, item in enumerate(ax.collections):
        if isinstance(item, matplotlib.collections.PolyCollection):
        # get path
            path, = item.get_paths()
            vertices = path.vertices

        # shift x-coordinates of path
            if not inner:
                if ii % 2: # -> to right
                    vertices[:,0] += delta
                else: # -> to left
                    vertices[:,0] -= delta
            else: # inner='box' adds another type of PollyCollection
                if ii % 3 == 0:
                    vertices[:,0] -= delta
                elif ii % 3 == 1:
                    vertices[:,0] += delta
                else: # ii % 3 = 2
                    pass

def violin_plot(metric):
    if metric != 'PEV':
        LR_braininf,mlrinf = snr_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected','inf');
        corr_LR_braininf,mcorr_lrinf = snr_brain(SNRinf[metric+'_corr_bx_x'],SNRinf[metric+'_sim_noise_x'],'Emp Corrected','inf');

        LR_brain100,mlr100 = snr_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','100');
        corr_LR_brain100,mcorr_lr100 = snr_brain(SNR100[metric+'_corr_bx_x'],SNR100[metric+'_sim_noise_x'],'Emp Corrected','100');

        LR_brain30,mlr30 = snr_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','30');
        corr_LR_brain30,mcorr_lr30 = snr_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','30');

        LR_brain10,mlr10 = snr_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','10');
        corr_LR_brain10,mcorr_lr10 = snr_brain(SNR10[metric+'_corr_bx_x'],SNR10[metric+'_sim_noise_x'],'Emp Corrected','10');

        aLR_braininf,mlrinf = snr_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected','ainf');
        acorr_LR_braininf,mcorr_lrinf = snr_brain(SNRinf[metric+'_corr_sm_fy_x'],SNRinf[metric+'_sim_noise_x'],'Emp Corrected','ainf');
        
        aLR_brain100,mlr100 = snr_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','a100');
        acorr_LR_brain100,mcorr_lr100 = snr_brain(SNR100[metric+'_corr_sm_fy_x'],SNR100[metric+'_sim_noise_x'],'Emp Corrected','a100');

        aLR_brain30,mlr30 = snr_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','a30');
        acorr_LR_brain30,mcorr_lr30 = snr_brain(SNR30[metric+'_corr_sm_fy_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','a30');
        
        aLR_brain10,mlr10 = snr_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','a10');
        acorr_LR_brain10,mcorr_lr10 = snr_brain(SNR10[metric+'_corr_sm_fy_x'],SNR10[metric+'_sim_noise_x'],'Emp Corrected','a10');
        pd_data = pd.concat([LR_braininf,corr_LR_braininf,LR_brain100,corr_LR_brain100,LR_brain30,corr_LR_brain30,LR_brain10,corr_LR_brain10,aLR_braininf,acorr_LR_braininf,aLR_brain100,acorr_LR_brain100,aLR_brain30,acorr_LR_brain30,aLR_brain10,acorr_LR_brain10])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
    else:
        LR_braininf,mlrinf = pev_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected','inf');
        corr_LR_braininf,mcorr_lrinf = pev_brain(SNRinf[metric+'_corr_bx_x'],SNRinf[metric+'_sim_noise_x'],'Emp Corrected','inf');

        LR_brain100,mlr100 = pev_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','100');
        corr_LR_brain100,mcorr_lr100 = pev_brain(SNR100[metric+'_corr_bx_x'],SNR100[metric+'_sim_noise_x'],'Emp Corrected','100');

        LR_brain30,mlr30 = pev_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','30');
        corr_LR_brain30,mcorr_lr30 = pev_brain(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','30');

        LR_brain10,mlr10 = pev_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','10');
        corr_LR_brain10,mcorr_lr10 = pev_brain(SNR10[metric+'_corr_bx_x'],SNR10[metric+'_sim_noise_x'],'Emp Corrected','10');

        aLR_braininf,mlrinf = pev_brain(SNRinf[metric+'_sim_corpt_x'],SNRinf[metric+'_true_x'],'Uncorrected','ainf');
        acorr_LR_braininf,mcorr_lrinf = pev_brain(SNRinf[metric+'_corr_sm_fy_x'],SNRinf[metric+'_sim_noise_x'],'Emp Corrected','ainf');

        aLR_brain100,mlr100 = pev_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','a100');
        acorr_LR_brain100,mcorr_lr100 = pev_brain(SNR100[metric+'_corr_sm_fy_x'],SNR100[metric+'_sim_noise_x'],'Emp Corrected','a100');

        aLR_brain30,mlr30 = pev_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','a30');
        acorr_LR_brain30,mcorr_lr30 = pev_brain(SNR30[metric+'_corr_sm_fy_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','a30');

        aLR_brain10,mlr10 = pev_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','a10');
        acorr_LR_brain10,mcorr_lr10 = pev_brain(SNR10[metric+'_corr_sm_fy_x'],SNR10[metric+'_sim_noise_x'],'Emp Corrected','a10');
        pd_data = pd.concat([LR_braininf,corr_LR_braininf,LR_brain100,corr_LR_brain100,LR_brain30,corr_LR_brain30,LR_brain10,corr_LR_brain10,aLR_braininf,acorr_LR_braininf,aLR_brain100,acorr_LR_brain100,aLR_brain30,acorr_LR_brain30,aLR_brain10,acorr_LR_brain10])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
    return pd_data

plt.subplots(3,1,figsize=(15,20))
sns.set(font_scale = 1.3)
sns.set_style("white")
palette = {'Uncorrected': 'crimson', 'Emp Corrected': 'limegreen'}
plt.subplot(3,1,1)
fa_pd_data = violin_plot('FA')
#ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=1500,palette=palette,inner="quart")
ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'H',y='value',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-10,100])
ax.legend(loc='upper right')
split_voilin()
ax.set_ylabel('APE in FA (%)',fontsize=15)
ax.set_xticklabels(['Inf','100','30','10','Inf','100','30','10'], fontsize=15)
plt.grid()

plt.subplot(3,1,2)
md_pd_data = violin_plot('MD')
ax = sns.violinplot(data=md_pd_data, hue = 'label', x = 'H',y='value',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-15,100])
ax.get_legend().remove()
split_voilin()
ax.set_ylabel('APE in MD (%)',fontsize=15)
ax.set_xticklabels(['Inf','100','30','10','Inf','100','30','10'], fontsize=15)
plt.grid()

plt.subplot(3,1,3)
v1_pd_data = violin_plot('PEV')
ax = sns.violinplot(data=v1_pd_data, hue = 'label', x = 'H',y='value',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-10,100])
ax.get_legend().remove()
split_voilin()
ax.set_ylabel('APE in PEV (%)',fontsize=15)
ax.set_xticklabels(['Inf','100','30','10','Inf','100','30','10'], fontsize=15)
plt.grid()

plt.tight_layout()
plt.show()
