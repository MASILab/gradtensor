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

SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNRinf.mat')
SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR30.mat')
SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR10.mat')

def snr_brain(FA_sim_noise_x,FA_true_x,label,x):
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

LR_braininf,mlrinf = snr_brain(SNRinf['FA_sim_corpt_x'],SNRinf['FA_true_x'],'Uncorrected','inf');
corr_LR_braininf,mcorr_lrinf = snr_brain(SNRinf['FA_corr_bx_x'],SNRinf['FA_sim_noise_x'],'Corrected','inf');

LR_brain100,mlr100 = snr_brain(SNR100['FA_sim_corpt_x'],SNR100['FA_true_x'],'Uncorrected','100');
corr_LR_brain100,mcorr_lr100 = snr_brain(SNR100['FA_corr_bx_x'],SNR100['FA_sim_noise_x'],'Corrected','100');

LR_brain30,mlr30 = snr_brain(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected','30');
corr_LR_brain30,mcorr_lr30 = snr_brain(SNR30['FA_corr_bx_x'],SNR30['FA_sim_noise_x'],'Corrected','30');

LR_brain10,mlr10 = snr_brain(SNR10['FA_sim_corpt_x'],SNR10['FA_true_x'],'Uncorrected','10');
corr_LR_brain10,mcorr_lr10 = snr_brain(SNR10['FA_corr_bx_x'],SNR10['FA_sim_noise_x'],'Corrected','10');

aLR_braininf,mlrinf = snr_brain(SNRinf['FA_sim_corpt_x'],SNRinf['FA_true_x'],'Uncorrected','ainf');
acorr_LR_braininf,mcorr_lrinf = snr_brain(SNRinf['FA_corr_sm_fy_x'],SNRinf['FA_sim_noise_x'],'Corrected','ainf');

aLR_brain100,mlr100 = snr_brain(SNR100['FA_sim_corpt_x'],SNR100['FA_true_x'],'Uncorrected','a100');
acorr_LR_brain100,mcorr_lr100 = snr_brain(SNR100['FA_corr_sm_fy_x'],SNR100['FA_sim_noise_x'],'Corrected','a100');

aLR_brain30,mlr30 = snr_brain(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected','a30');
acorr_LR_brain30,mcorr_lr30 = snr_brain(SNR30['FA_corr_sm_fy_x'],SNR30['FA_sim_noise_x'],'Corrected','a30');

aLR_brain10,mlr10 = snr_brain(SNR10['FA_sim_corpt_x'],SNR10['FA_true_x'],'Uncorrected','a10');
acorr_LR_brain10,mcorr_lr10 = snr_brain(SNR10['FA_corr_sm_fy_x'],SNR10['FA_sim_noise_x'],'Corrected','a10');

fig, ax = plt.subplots(1,1)
sns.set(font_scale = 2)
pd_data = pd.concat([LR_braininf,corr_LR_braininf,LR_brain100,corr_LR_brain100,LR_brain30,corr_LR_brain30,LR_brain10,corr_LR_brain10,aLR_braininf,acorr_LR_braininf,aLR_brain100,acorr_LR_brain100,aLR_brain30,acorr_LR_brain30,aLR_brain10,acorr_LR_brain10])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
#pd_data = pd_data.melt(id_vars=['label'], value_vars=['inf','100','30','10'],var_name='H', value_name='value')
ax = sns.violinplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',gridsize=3000,split=True,dodge=False,inner=None,linewidth=2,scale="width")
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
ax.set_xticklabels(['Inf','100','30','10','Inf','100','30','10'])
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(loc='upper left')
plt.title('Comparison of correction technique on FA in tensor simulation (GM)')
#plt.ylim([-0.1,2])
plt.tight_layout()
plt.grid()
plt.show()
