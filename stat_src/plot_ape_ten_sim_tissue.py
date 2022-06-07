import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.collections
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy.io
import h5py
import mat73

#SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNRinf.mat')
#SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR30.mat')
#SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_50_30d_SNR10.mat')

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

def snr_csf(FA_sim_noise_x,FA_true_x,label,x,pev=False):
    if pev:
        err_fa = angular_error(FA_sim_noise_x,FA_true_x)
        pe_fa = err_fa
    else:
        err_fa = FA_sim_noise_x - FA_true_x
        pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(pe_fa)
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[0,]
    df_brain = pd.DataFrame(brain).assign(label = label)
    df_brain = df_brain.assign(x = x)
    #idf_brain = df_brain.rename(columns={0:x})
    median = df_brain.groupby(['label'])[0].median().values
    return df_brain, median

def snr_wm(FA_sim_noise_x,FA_true_x,label,x,pev=False):
    if pev:
        err_fa = angular_error(FA_sim_noise_x,FA_true_x)
        pe_fa = err_fa
    else:
        err_fa = FA_sim_noise_x - FA_true_x
        pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(pe_fa)
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[7,]
    df_brain = pd.DataFrame(brain).assign(label = label)
    df_brain = df_brain.assign(x = x)
    #idf_brain = df_brain.rename(columns={0:x})
    median = df_brain.groupby(['label'])[0].median().values
    return df_brain, median

def snr_gm(FA_sim_noise_x,FA_true_x,label,x,pev=False):
    if pev:
        err_fa = angular_error(FA_sim_noise_x,FA_true_x)
        pe_fa = err_fa
    else:
        err_fa = FA_sim_noise_x - FA_true_x
        pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(pe_fa)
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[2,]
    df_brain = pd.DataFrame(brain).assign(label = label)
    df_brain = df_brain.assign(x = x)
    #idf_brain = df_brain.rename(columns={0:x})
    median = df_brain.groupby(['label'])[0].median().values
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
        LR_brain30,mlr30 = snr_csf(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','CSF');
        corr_LR_brain30,mcorr_lr30 = snr_csf(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','CSF');

        gLR_brain30,mlr30 = snr_gm(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','GM');
        gcorr_LR_brain30,mcorr_lr30 = snr_gm(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','GM');

        wLR_brain30,mlr30 = snr_wm(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','WM');
        wcorr_LR_brain30,mcorr_lr30 = snr_wm(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','WM');
        pd_data = pd.concat([wLR_brain30,wcorr_LR_brain30,gLR_brain30,gcorr_LR_brain30,LR_brain30,corr_LR_brain30])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
    else:
        LR_brain30,mlr30 = snr_csf(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','CSF',True);
        corr_LR_brain30,mcorr_lr30 = snr_csf(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','CSF',True);

        gLR_brain30,mlr30 = snr_gm(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','GM',True);
        gcorr_LR_brain30,mcorr_lr30 = snr_gm(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','GM',True);

        wLR_brain30,mlr30 = snr_wm(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','WM',True);
        wcorr_LR_brain30,mcorr_lr30 = snr_wm(SNR30[metric+'_corr_bx_x'],SNR30[metric+'_sim_noise_x'],'Emp Corrected','WM',True);
        pd_data = pd.concat([wLR_brain30,wcorr_LR_brain30,gLR_brain30,gcorr_LR_brain30,LR_brain30,corr_LR_brain30])
        pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
    return pd_data

plt.subplots(3,1,figsize=(15,20))
sns.set(font_scale = 1.3)
sns.set_style("white")
palette = {'Uncorrected': 'crimson', 'Emp Corrected': 'limegreen'}
plt.subplot(3,1,1)
fa_pd_data = violin_plot('FA')
#ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=1500,palette=palette,inner="quart")
ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-10,100])
ax.legend(loc='upper right')
split_voilin()
ax.set_ylabel('APE in FA (%)',fontsize=15)
ax.set_xticklabels(['WM','GM','CSF'],fontsize=15)
plt.grid()

plt.subplot(3,1,2)
md_pd_data = violin_plot('MD')
ax = sns.violinplot(data=md_pd_data, hue = 'label',  x = 'x',y='Abs percent error (%)',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-15,100])
ax.get_legend().remove()
split_voilin()
ax.set_ylabel('APE in MD (%)',fontsize=15)
ax.set_xticklabels(['WM','GM','CSF'],fontsize=15)
plt.grid()

plt.subplot(3,1,3)
v1_pd_data = violin_plot('PEV')
ax = sns.violinplot(data=v1_pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',gridsize=4000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
ax.set_ylim([-10,100])
ax.get_legend().remove()
split_voilin()
ax.set_ylabel('APE in PEV (%)',fontsize=15)
ax.set_xticklabels(['WM','GM','CSF'],fontsize=15)
plt.grid()

plt.tight_layout()
plt.show()
