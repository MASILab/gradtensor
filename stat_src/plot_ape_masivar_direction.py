import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import pandas as pd
import math


mask32 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/mask.nii').get_fdata()
mask40 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/10/mask.nii').get_fdata()
mask96 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/11/mask.nii').get_fdata()


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

def plot_one_snr(corpt_fa,true_fa,label,x,pev=False,mask=None):
    if pev:
        err_fa = angular_error(corpt_fa,true_fa)
        pe_fa = err_fa
    else:
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
    true_fa32 =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+metric+'.nii').get_fdata()
    noise_fa32 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa32 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa32 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa32 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_1/approx_corrected_'+metric+'.nii').get_fdata()

    true_fa40 =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/10/true_'+metric+'.nii').get_fdata()
    noise_fa40 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d40_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa40 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d40_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa40 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d40_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa40 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d40_1/approx_corrected_'+metric+'.nii').get_fdata()

    true_fa96 =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/11/true_'+metric+'.nii').get_fdata()
    noise_fa96 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d96_1/uncorrected_Nest_'+metric+'.nii').get_fdata()
    corpt_fa96 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d96_1/uncorrected_'+metric+'.nii').get_fdata()
    bx_corr_fa96 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d96_1/emp/emp_corrected_'+metric+'.nii').get_fdata()
    sm_corr_fa96 = nib.load('/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d96_1/approx_corrected_'+metric+'.nii').get_fdata()

    if metric != 'primary_eigvec':
        lr_noise30 = plot_one_snr(corpt_fa32,true_fa32,'Uncorrected','d32',mask=mask32)
        corr_lr30 = plot_one_snr(bx_corr_fa32,noise_fa32,'Emp Corrected','d32',mask=mask32)

        lr_noise30_40 = plot_one_snr(corpt_fa40,true_fa40,'Uncorrected','d40',mask=mask40)
        corr_lr30_40 = plot_one_snr(bx_corr_fa40,noise_fa40,'Emp Corrected','d40',mask=mask40)

        lr_noise30_96 = plot_one_snr(corpt_fa96,true_fa96,'Uncorrected','d96',mask=mask96)
        corr_lr30_96 = plot_one_snr(bx_corr_fa96,noise_fa96,'Emp Corrected','d96',mask=mask96)
        pd_data = pd.concat([lr_noise30_96,corr_lr30_96,lr_noise30_40,corr_lr30_40,lr_noise30,corr_lr30])
        pd_data = pd_data.melt(id_vars=['label'], value_vars=['d32','d40','d96'],var_name='H', value_name='value')
    else:
        lr_noise30 = plot_one_snr(corpt_fa32,true_fa32,'Uncorrected','d32',True,mask=mask32)
        corr_lr30 = plot_one_snr(bx_corr_fa32,noise_fa32,'Emp Corrected','d32',True,mask=mask32)

        lr_noise30_40 = plot_one_snr(corpt_fa40,true_fa40,'Uncorrected','d40',True,mask=mask40)
        corr_lr30_40 = plot_one_snr(bx_corr_fa40,noise_fa40,'Emp Corrected','d40',True,mask=mask40)

        lr_noise30_96 = plot_one_snr(corpt_fa96,true_fa96,'Uncorrected','d96',True,mask=mask96)
        corr_lr30_96 = plot_one_snr(bx_corr_fa96,noise_fa96,'Emp Corrected','d96',True,mask=mask96)
        pd_data = pd.concat([lr_noise30_96,corr_lr30_96,lr_noise30_40,corr_lr30_40,lr_noise30,corr_lr30])
        pd_data = pd_data.melt(id_vars=['label'], value_vars=['d32','d40','d96'],var_name='H', value_name='value')
    return pd_data

plt.subplots(3,1,figsize=(15,20))
sns.set(font_scale = 1.3)
sns.set_style("white")
palette = {'Uncorrected': 'crimson', 'Emp Corrected': 'limegreen'}
plt.subplot(3,1,1)
fa_pd_data = violin_plot('fa')
ax = sns.violinplot(data=fa_pd_data, hue = 'label', x = 'H',y='value',gridsize=5000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette,bw=0.065)
#ax.set_ylim([-1,4])
ax.set_ylim([-10,100])
ax.legend(loc='upper right')
ax.set_xlabel(' ')
split_voilin()
#ax.set(xlabel=' ',ylabel='APE in FA (%)',fontsize= 22)
ax.set_ylabel('APE in FA (%)',fontsize=15)
ax.set_xticklabels(['32','40','96'],fontsize=15)
plt.grid()

plt.subplot(3,1,2)
md_pd_data = violin_plot('md')
ax = sns.violinplot(data=md_pd_data, hue = 'label', x = 'H',y='value',gridsize=10000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette,bw=0.06)
#ax.set_ylim([-1,4])
ax.set_ylim([-10,100])
ax.get_legend().remove()
split_voilin()
ax.set_ylabel('APE in MD (%)',fontsize=15)
ax.set_xlabel(' ')
ax.set_xticklabels(['32','40','96'],fontsize=15)
plt.grid()

plt.subplot(3,1,3)
v1_pd_data = violin_plot('primary_eigvec')
ax = sns.violinplot(data=v1_pd_data, hue = 'label', x = 'H',y='value',gridsize=5000,split=True,dodge=False,linewidth=2,scale="width",inner=None,palette=palette)
#ax.set_ylim([-1,4])
ax.set_ylim([-10,100])
ax.get_legend().remove()
split_voilin()
ax.set_xlabel(' ')
ax.set_ylabel('APE in PEV (%)',fontsize=15)
ax.set_xticklabels(['32','40','96'],fontsize=15)
plt.grid()  

plt.tight_layout()
plt.show()
