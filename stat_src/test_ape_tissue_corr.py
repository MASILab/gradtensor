import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import pandas as pd
import math

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__fa.nii').get_fdata()
noise_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
corpt_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_fa.nii').get_fdata()
bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/approx_corrected_fa.nii').get_fdata()
seg =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/seg_reg_skull_stripped.nii.gz').get_fdata()

def plot_csf(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(seg == 1,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    df_ape_fa = df_ape_fa.assign(x = x)
    #df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa

def plot_gm(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(seg == 2,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    df_ape_fa = df_ape_fa.assign(x = x)
    #df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa

def plot_wm(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(seg == 3,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    df_ape_fa = df_ape_fa.assign(x = x)
    #df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa

lr_noise30 = plot_csf(corpt_fa30,true_fa,'Uncorrected','CSF')
corr_lr30 = plot_csf(bx_corr_fa30,noise_fa30,'Corrected','CSF')

glr_noise30 = plot_gm(corpt_fa30,true_fa,'Uncorrected','GM')
gcorr_lr30 = plot_gm(bx_corr_fa30,noise_fa30,'Corrected','GM')

wlr_noise30 = plot_wm(corpt_fa30,true_fa,'Uncorrected','WM')
wcorr_lr30 = plot_wm(bx_corr_fa30,noise_fa30,'Corrected','WM')

fig, ax = plt.subplots(1,1)
sns.set(font_scale = 2)
pd_data = pd.concat([wlr_noise30,wcorr_lr30,glr_noise30,gcorr_lr30,lr_noise30,corr_lr30])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
ax = sns.violinplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',gridsize=2000,split=True,dodge=False,inner=None,linewidth=2,scale="width")
inner=None
# offset stuff
delta = 0.02
for ii, item in enumerate(ax.collections):
    # axis contains PolyCollections and PathCollections
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
#sns.despine(left=True)
ax.set_xlabel('SNR',fontsize = 17)
ax.set_ylabel('Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )',fontsize = 17)
#ax.set_xticklabels(['WM','GM''CSF'])
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(loc='upper left')
plt.title('Corruption and Emperical Correction on FA in synethic data for different tissue at SNR 30')
plt.ylim([-0.5,4])
plt.tight_layout()
plt.grid()
plt.show()
