import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
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
    pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    #df_ape_fa = df_ape_fa.assign(x = x)
    df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa

def plot_one_snr_pev(corpt_fa,true_fa,label,x):
    err_fa = angular_error(corpt_fa,true_fa)
    #pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(err_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    #df_ape_fa = df_ape_fa.assign(x = x)
    df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa


def split_voilin(delta=0.03,inner=None):
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
    #true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__'+metric+'.nii').get_fdata()
    true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+metric+'.nii').get_fdata()
    noise_fainf = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fainf = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fainf = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fainf = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1/approx_corrected_'+metric+'.nii').get_fdata()

    noise_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR100_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR100_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR100_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR100_d32_1/approx_corrected_'+metric+'.nii').get_fdata()

    noise_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_corrected_'+metric+'.nii').get_fdata()

    noise_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR10_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR10_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR10_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR10_d32_1/approx_corrected_'+metric+'.nii').get_fdata()
    if metric != 'primary_eigvec':
        lr_noiseinf = plot_one_snr(corpt_fainf,true_fa,'Uncorrected','emp_inf')
        corr_lrinf = plot_one_snr(bx_corr_fainf,noise_fainf,'Emp Corrected','emp_inf')

        lr_noise100 = plot_one_snr(corpt_fa100,true_fa,'Uncorrected','emp_100')
        corr_lr100 = plot_one_snr(bx_corr_fa100,noise_fa100,'Emp Corrected','emp_100')
    
        lr_noise30 = plot_one_snr(corpt_fa30,true_fa,'Uncorrected','emp_30')
        corr_lr30 = plot_one_snr(bx_corr_fa30,noise_fa30,'Emp Corrected','emp_30')

        lr_noise10 = plot_one_snr(corpt_fa10,true_fa,'Uncorrected','emp_10')
        corr_lr10 = plot_one_snr(bx_corr_fa10,noise_fa10,'Emp Corrected','emp_10')

        alr_noiseinf = plot_one_snr(corpt_fainf,true_fa,'Uncorrected','app_inf')
        acorr_lrinf = plot_one_snr(sm_corr_fainf,noise_fainf,'Approx Corrected','app_inf')

        alr_noise100 = plot_one_snr(corpt_fa100,true_fa,'Uncorrected','app_100')
        acorr_lr100 = plot_one_snr(sm_corr_fa100,noise_fa100,'Approx Corrected','app_100')

        alr_noise30 = plot_one_snr(corpt_fa30,true_fa,'Uncorrected','app_30')
        acorr_lr30 = plot_one_snr(sm_corr_fa30,noise_fa30,'Approx Corrected','app_30')

        alr_noise10 = plot_one_snr(corpt_fa10,true_fa,'Uncorrected','app_10')
        acorr_lr10 = plot_one_snr(sm_corr_fa10,noise_fa10,'Approx Corrected','app_10')

    else:
        lr_noiseinf = plot_one_snr_pev(corpt_fainf,true_fa,'Uncorrected','emp_inf')
        corr_lrinf = plot_one_snr_pev(bx_corr_fainf,noise_fainf,'Emp Corrected','emp_inf')

        lr_noise100 = plot_one_snr_pev(corpt_fa100,true_fa,'Uncorrected','emp_100')
        corr_lr100 = plot_one_snr_pev(bx_corr_fa100,noise_fa100,'Emp Corrected','emp_100')

        lr_noise30 = plot_one_snr_pev(corpt_fa30,true_fa,'Uncorrected','emp_30')
        corr_lr30 = plot_one_snr_pev(bx_corr_fa30,noise_fa30,'Emp Corrected','emp_30')

        lr_noise10 = plot_one_snr_pev(corpt_fa10,true_fa,'Uncorrected','emp_10')
        corr_lr10 = plot_one_snr_pev(bx_corr_fa10,noise_fa10,'Emp Corrected','emp_10')

        alr_noiseinf = plot_one_snr_pev(corpt_fainf,true_fa,'Uncorrected','app_inf')
        acorr_lrinf = plot_one_snr_pev(sm_corr_fainf,noise_fainf,'Approx Corrected','app_inf')

        alr_noise100 = plot_one_snr_pev(corpt_fa100,true_fa,'Uncorrected','app_100')
        acorr_lr100 = plot_one_snr_pev(sm_corr_fa100,noise_fa100,'Approx Corrected','app_100')

        alr_noise30 = plot_one_snr_pev(corpt_fa30,true_fa,'Uncorrected','app_30')
        acorr_lr30 = plot_one_snr_pev(sm_corr_fa30,noise_fa30,'Approx Corrected','app_30')

        alr_noise10 = plot_one_snr_pev(corpt_fa10,true_fa,'Uncorrected','app_10')
        acorr_lr10 = plot_one_snr_pev(sm_corr_fa10,noise_fa10,'Approx Corrected','app_10')
    return lr_noiseinf, lr_noise100, lr_noise30, lr_noise10, corr_lrinf, corr_lr100, corr_lr30, corr_lr10, alr_noiseinf, alr_noise100, alr_noise30, alr_noise10, acorr_lrinf, acorr_lr100, acorr_lr30, acorr_lr10

lr_noiseinf, lr_noise100, lr_noise30, lr_noise10, corr_lrinf, corr_lr100, corr_lr30, corr_lr10, alr_noiseinf, alr_noise100, alr_noise30, alr_noise10, acorr_lrinf, acorr_lr100, acorr_lr30, acorr_lr10 = violin_plot('fa')
def cohens_d(c1,c0):
    cohens_d = (mean(c0) - mean(c1)) / (sqrt((stdev(c0) ** 2 + stdev(c1) ** 2) / 2))
    return cohens_d

def run_cohens_d(m):
    lr_noiseinf, lr_noise100, lr_noise30, lr_noise10, corr_lrinf, corr_lr100, corr_lr30, corr_lr10, alr_noiseinf, alr_noise100, alr_noise30, alr_noise10, acorr_lrinf, acorr_lr100, acorr_lr30, acorr_lr10 = violin_plot(m)

    corpt_emp_inf = cohens_d(lr_noiseinf['emp_inf'],corr_lrinf['emp_inf'])
    corpt_emp_100 = cohens_d(lr_noise100['emp_100'],corr_lr100['emp_100'])
    corpt_emp_30 = cohens_dlr_noise30['emp_30'],corr_lr30['emp_30'])
    corpt_emp_10 = cohens_d(lr_noise10['emp_10'],corr_lr10['emp_10'])

    corpt_app_inf = cohens_d(lr_noiseinf['emp_inf'],acorr_lrinf['app_inf'])
    corpt_app_100 = cohens_d(lr_noise100['emp_100'],acorr_lr100['app_100'])
    corpt_app_30 = cohens_d(lr_noise30['emp_30'],acorr_lr30['app_30'])
    corpt_app_10 = cohens_d(lr_noise10['emp_10'],acorr_lr10['app_10'])

    emp_app_inf = cohens_d(corr_lrinf['emp_inf'],acorr_lrinf['app_inf'])
    emp_app_100 = cohens_d(corr_lr100['emp_100'],acorr_lr100['app_100'])
    emp_app_30 = cohens_d(corr_lr30['emp_30'],acorr_lr30['app_30'])
    emp_app_10 = cohens_d(corr_lr10['emp_10'],acorr_lr10['app_10'])
    print(corpt_emp_inf,corpt_emp_100,corpt_emp_30,corpt_emp_10,corpt_app_inf,corpt_app_100,corpt_app_30,corpt_app_10,emp_app_inf,emp_app_100,emp_app_30,emp_app_10)

run_cohens_d(sys.args[1])
