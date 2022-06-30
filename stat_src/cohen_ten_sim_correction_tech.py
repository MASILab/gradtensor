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
import sys

SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNRinf.mat')
SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR30.mat')
SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR10.mat')

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
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[1,] # WM13
    df_brain = pd.DataFrame(brain).assign(x = x)
    df_brain = df_brain.assign(label = label)
    median = df_brain.groupby(['label', 'x'])[0].median().values
    return df_brain, median

def pev_brain(FA_sim_noise_x,FA_true_x,label,x):
    err_fa = angular_error(FA_sim_noise_x,FA_true_x)
    #pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(err_fa)
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[1,] # WM
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
        acorr_LR_braininf,mcorr_lrinf = snr_brain(SNRinf[metric+'_corr_sm_fy_x'],SNRinf[metric+'_sim_noise_x'],'Approx Corrected','ainf');
        
        aLR_brain100,mlr100 = snr_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','a100');
        acorr_LR_brain100,mcorr_lr100 = snr_brain(SNR100[metric+'_corr_sm_fy_x'],SNR100[metric+'_sim_noise_x'],'Approx Corrected','a100');

        aLR_brain30,mlr30 = snr_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','a30');
        acorr_LR_brain30,mcorr_lr30 = snr_brain(SNR30[metric+'_corr_sm_fy_x'],SNR30[metric+'_sim_noise_x'],'Approx Corrected','a30');
        
        aLR_brain10,mlr10 = snr_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','a10');
        acorr_LR_brain10,mcorr_lr10 = snr_brain(SNR10[metric+'_corr_sm_fy_x'],SNR10[metric+'_sim_noise_x'],'Approx Corrected','a10');

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
        acorr_LR_braininf,mcorr_lrinf = pev_brain(SNRinf[metric+'_corr_sm_fy_x'],SNRinf[metric+'_sim_noise_x'],'Approx Corrected','ainf');

        aLR_brain100,mlr100 = pev_brain(SNR100[metric+'_sim_corpt_x'],SNR100[metric+'_true_x'],'Uncorrected','a100');
        acorr_LR_brain100,mcorr_lr100 = pev_brain(SNR100[metric+'_corr_sm_fy_x'],SNR100[metric+'_sim_noise_x'],'Approx Corrected','a100');

        aLR_brain30,mlr30 = pev_brain(SNR30[metric+'_sim_corpt_x'],SNR30[metric+'_true_x'],'Uncorrected','a30');
        acorr_LR_brain30,mcorr_lr30 = pev_brain(SNR30[metric+'_corr_sm_fy_x'],SNR30[metric+'_sim_noise_x'],'Approx Corrected','a30');

        aLR_brain10,mlr10 = pev_brain(SNR10[metric+'_sim_corpt_x'],SNR10[metric+'_true_x'],'Uncorrected','a10');
        acorr_LR_brain10,mcorr_lr10 = pev_brain(SNR10[metric+'_corr_sm_fy_x'],SNR10[metric+'_sim_noise_x'],'Approx Corrected','a10');

    return LR_braininf,LR_brain100,LR_brain30,LR_brain10,corr_LR_braininf,corr_LR_brain100,corr_LR_brain30,corr_LR_brain10,aLR_braininf,aLR_brain100,aLR_brain30,aLR_brain10

def cohens_d(c1,c0):
    cohens_d = (mean(c0) - mean(c1)) / (sqrt((stdev(c0) ** 2 + stdev(c1) ** 2) / 2))
    return cohens_d

def run_cohens_d(m):
    LR_braininf,LR_brain100,LR_brain30,LR_brain10,corr_LR_braininf,corr_LR_brain100,corr_LR_brain30,corr_LR_brain10,aLR_braininf,aLR_brain100,aLR_brain30,aLR_brain10 = violin_plot(m)
    corpt_emp_inf = cohens_d(LR_braininf['inf'],corr_LR_braininf['inf'])
    corpt_emp_100 = cohens_d(LR_brain100['100'],corr_LR_brain100['100'])
    corpt_emp_30 = cohens_d(LR_brain30['30'],corr_LR_brain30['30'])
    corpt_emp_10 = cohens_d(LR_brain10['10'],corr_LR_brain10['10'])

    corpt_app_inf = cohens_d(LR_braininf['inf'],acorr_LR_braininf['ainf'])
    corpt_app_100 = cohens_d(LR_brain100['100'],acorr_LR_brain100['a100'])
    corpt_app_30 = cohens_d(LR_brain30['30'],acorr_LR_brain30['a30'])
    corpt_app_10 = cohens_d(LR_brain10['10'],acorr_LR_brain10['a10'])

    emp_app_inf = cohens_d(corr_LR_braininf['inf'],acorr_LR_braininf['ainf'])
    emp_app_100 = cohens_d(corr_LR_brain100['100'],acorr_LR_brain100['a100'])
    emp_app_30 = cohens_d(corr_LR_brain30['30'],acorr_LR_brain30['a30'])
    emp_app_10 = cohens_d(corr_LR_brain10['10'],acorr_LR_brain10['a10'])
    print(corpt_emp_inf,corpt_emp_100,corpt_emp_30,corpt_emp_10,corpt_app_inf,corpt_app_100,corpt_app_30,corpt_app_10,emp_app_inf,emp_app_100,emp_app_30,emp_app_10)

run_cohens_d(sys.argv[1])
