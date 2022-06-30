import scipy.io
import h5py
import mat73
import pandas as pd
import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math

SNRinf = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNRinf.mat')
#SNR100 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR100.mat')
SNR30 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR30.mat')
#SNR10 = mat73.loadmat('/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_10_30d_SNR10.mat')

def snr_brain(FA_sim_noise_x,FA_true_x,label,x):
    err_fa = FA_sim_noise_x - FA_true_x
    pe_fa =  (err_fa / FA_true_x) * 100
    ape_fa = np.abs(pe_fa)
    cube = np.reshape(ape_fa, [10, 10, 19406])
    zzz = np.squeeze(np.mean(cube,2))
    zzz_transp = np.transpose(zzz)
    brain = zzz_transp[2,]
    df_brain = pd.DataFrame(brain).assign(x = x)
    df_brain = df_brain.assign(label = label)
    median = df_brain.groupby(['label', 'x'])[0].median().values
    return df_brain, median

noise_braininf,minf_noise= snr_brain(SNRinf['FA_sim_noise_x'],SNRinf['FA_true_x'],'Uncorrected SNR=inf','inf_noise');
noiseLR_braininf,minf_lr_noise = snr_brain(SNRinf['FA_sim_corpt_x'],SNRinf['FA_true_x'],'Uncorrected SNR=inf','inf_lr_noise');
LR_braininf,minf_lr = snr_brain(SNRinf['FA_sim_corpt_x'],SNRinf['FA_sim_noise_x'],'Uncorrected SNR=inf','inf_lr');

noise_brain30,mnoise = snr_brain(SNR30['FA_sim_noise_x'],SNR30['FA_true_x'],'Uncorrected SNR=30','noise');
noiseLR_brain30,mlr_noise = snr_brain(SNR30['FA_sim_corpt_x'],SNR30['FA_true_x'],'Uncorrected SNR=30','lr_noise');
LR_brain30,mlr = snr_brain(SNRinf['FA_sim_corpt_x'],SNR30['FA_sim_noise_x'],'Uncorrected SNR=30','lr');

corr_noise_brain30,mcorr_noise = snr_brain(SNR30['FA_corr_noise_bx_x'],SNR30['FA_true_x'],'Emp corrected SNR=30','corr_noise');
corr_noiseLR_brain30,mcorr_lr_noise = snr_brain(SNR30['FA_corr_bx_x'],SNR30['FA_true_x'],'Emp corrected SNR=30','corr_lr_noise');
corr_LR_brain30,mcorr_lr = snr_brain(SNR30['FA_corr_bx_x'],SNR30['FA_sim_noise_x'],'Emp corrected SNR=30','corr_lr');

#plt.figure()
sns.set(font_scale = 1.3)
pd_data = pd.concat([noise_braininf,noise_brain30,corr_noise_brain30,noiseLR_braininf,noiseLR_brain30,corr_noiseLR_brain30,LR_braininf,LR_brain30,corr_LR_brain30])
pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
sns.set_style("white")
ax = sns.violinplot(data=pd_data, hue = 'label', x = 'x',y='Abs percent error (%)',dodge=False,width=.5,gridsize=2000,inner="stick")
#sns.despine(left=True)

ax.set(xlabel=' ', ylabel = 'Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )')

ax.set_xticklabels(['','noise FA - GT FA','','','noise+LR FA - GT FA','','','noise+LR FA - noise FA',''])

m1 = np.concatenate([minf_noise, mnoise, mcorr_noise, minf_lr_noise, mlr_noise,mcorr_lr_noise,minf_lr, mlr, mcorr_lr])
print(m1)
mL1 = [(np.round(s, 2)) for s in m1]
print(mL1)
for xtick in ax.get_xticks():
    ax.text(xtick-.2,mL1[xtick] + ( mL1[xtick]*0.05),mL1[xtick],
            horizontalalignment='center',size='x-small',color='k',weight='semibold')

plt.legend(loc='upper right')
plt.title('Effect of noise and LR on FA in tensor simulation (CSF)')
plt.ylim([-15,150])
plt.grid()
plt.show()



