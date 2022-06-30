import scipy.io
import h5py
import mat73
import pandas as pd
import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import math

#SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNRinf.mat')
#SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR30.mat')
#SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR10.mat')

def snr_csf(FA_sim_noise_x,FA_true_x,label,x):
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

def snr_wm(FA_sim_noise_x,FA_true_x,label,x):
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

def snr_gm(FA_sim_noise_x,FA_true_x,label,x):
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


LR_brain30,mlr30 = snr_csf(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected','CSF');
corr_LR_brain30,mcorr_lr30 = snr_csf(SNR30['FA_corr_bx_x'],SNR30['FA_sim_noise_x'],'Corrected','CSF');

gLR_brain30,mlr30 = snr_gm(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected','GM');
gcorr_LR_brain30,mcorr_lr30 = snr_gm(SNR30['FA_corr_bx_x'],SNR30['FA_sim_noise_x'],'Corrected','GM');

wLR_brain30,mlr30 = snr_wm(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected','WM');
wcorr_LR_brain30,mcorr_lr30 = snr_wm(SNR30['FA_corr_bx_x'],SNR30['FA_sim_noise_x'],'Corrected','WM');


fig, ax = plt.subplots(1,1)
sns.set(font_scale = 2)
pd_data = pd.concat([wLR_brain30,wcorr_LR_brain30,gLR_brain30,gcorr_LR_brain30,LR_brain30,corr_LR_brain30])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
#pd_data = pd_data.melt(id_vars=['label'], value_vars=['inf','100','30','10'],var_name='H', value_name='value')
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
plt.legend(loc='center right')
plt.title('Corruption and Emperical Correction on FA in tensor simulation for different tissue at SNR 30')
plt.ylim([-0.1,2])
plt.tight_layout()
plt.grid()
plt.show()
