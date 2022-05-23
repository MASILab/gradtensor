import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

true_fa_32 =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__md.nii').get_fdata()
corpt_fa100_32 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_md.nii').get_fdata()
sm_corr_fa100_32 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa100_32 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/emp_corrected_md.nii').get_fdata()
seg_32 =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/seg_reg_skull_stripped.nii.gz').get_fdata()

true_fa_40 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_/true__md.nii').get_fdata() 
corpt_fa100_40 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_1/uncorrected_md.nii').get_fdata()
sm_corr_fa100_40 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa100_40 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_1/emp_corrected_md.nii').get_fdata()
seg_40 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_1/seg_reg_skull_stripped.nii.gz').get_fdata()

true_fa_91 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d94_/true__md.nii').get_fdata()
corpt_fa100_91 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d94_1/uncorrected_md.nii').get_fdata()
sm_corr_fa100_91 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d94_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa100_91 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d94_1/emp_corrected_md.nii').get_fdata()
seg_91 =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d94_1/seg_reg_skull_stripped.nii.gz').get_fdata()

# wm
def plot_one_snr(corpt_fa,true_fa,acq,seg):
    err_fa = corpt_fa - true_fa
    pe_fa = 100 * err_fa / true_fa
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa_csf = np.where(seg == 1, ape_fa , math.nan) # for csf 
    ape_fa_wm = np.where(seg == 3, ape_fa , math.nan) 
    ape_fa_gm = np.where(seg == 2, ape_fa , math.nan)
    ape_fa_csf = ape_fa_csf.reshape(-1)
    ape_fa_gm = ape_fa_gm.reshape(-1)
    ape_fa_wm = ape_fa_wm.reshape(-1)
    ape_fa_csf = [x for x in ape_fa_csf if math.isnan(x) == False]
    ape_fa_gm = [x for x in ape_fa_gm if math.isnan(x) == False]
    ape_fa_wm = [x for x in ape_fa_wm if math.isnan(x) == False]
    df1 = pd.DataFrame(ape_fa_gm).assign(Tissue='GM')
    df2 = pd.DataFrame(ape_fa_wm).assign(Tissue='WM')
    df3 = pd.DataFrame(ape_fa_csf).assign(Tissue='CSF')
    cdf = pd.concat([df3, df1, df2])
    cdf = cdf.assign(ACQ = acq)
    median = cdf.groupby(['ACQ', 'Tissue'])[0].median().values
    idx = [0,1,2]
    median_sort = median[idx]
    return cdf, median

cdf32,median32 = plot_one_snr(corpt_fa100_32,true_fa_32,'32',seg_32)
cdf_ec32,median_ec32 = plot_one_snr(bx_corr_fa100_32,true_fa_32,'32 ec',seg_32)
cdf_ac32,median_ac32 = plot_one_snr(sm_corr_fa100_32,true_fa_32,'32 ac',seg_32)

cdf40,median40 = plot_one_snr(corpt_fa100_40,true_fa_40,'40',seg_40)
cdf_ec40,median_ec40 = plot_one_snr(bx_corr_fa100_40,true_fa_40,'40 ec',seg_40)
cdf_ac40,median_ac40 = plot_one_snr(sm_corr_fa100_40,true_fa_40,'40 ac',seg_40)

cdf91,median91 = plot_one_snr(corpt_fa100_91,true_fa_91,'91',seg_91)
cdf_ec91,median_ec91 = plot_one_snr(bx_corr_fa100_91,true_fa_91,'91 ec',seg_91)
cdf_ac91,median_ac91 = plot_one_snr(sm_corr_fa100_91,true_fa_91,'91 ac',seg_91)

all_cdf = pd.concat([cdf32,cdf_ec32,cdf_ac32,cdf40,cdf_ec40,cdf_ac40, cdf91, cdf_ec91,cdf_ac91])
all_cdf = all_cdf.rename(columns={0:'Abs percent error (%)'})
plt.figure()
sns.set(font_scale = 1.3)
flierprops = dict(markerfacecolor='0.75', markersize=2,
              linestyle='none')
testPlot = sns.boxplot(data=all_cdf,hue='Tissue',x = 'ACQ',y='Abs percent error (%)',palette="Set2", flierprops=flierprops)

#m1 = all_cdf.groupby(['SNR', 'Tissue'])['Abs percent error (%)'].median().values
m1 = np.concatenate([median32,median_ec32,median_ac32,median40,median_ec40,median_ac40,median91,median_ec91,median_ac91])
mL1 = [str(np.round(s, 2)) for s in m1]
ind = 0
for tick in range(len(testPlot.get_xticklabels())):
    testPlot.text(tick+.4, m1[ind+2]+1, mL1[ind+2],  horizontalalignment='center',  color='k', weight='semibold', fontsize=10)
    testPlot.text(tick+.1, m1[ind+1]+1, mL1[ind+1],  horizontalalignment='center',  color='k', weight='semibold', fontsize=10)
    testPlot.text(tick-.3, m1[ind]+1, mL1[ind], horizontalalignment='center', color='k', weight='semibold', fontsize=10)
    ind += 3 

#testPlot.set_xticklabels(['','32','','','40','','','91',''])
tick_labels=['Corruption','\n\n ','Emp Corr','\n\n32','Appr Corr','\n\n','Corruption','\n\n ','Emp Corr','\n\n40','Appr Corr','\n\n','Corruption','\n\n ','Emp Corr','\n\n91','Appr Corr','\n\n',]#'Corruption','\n\n ','Emp Corr','\n\n10','Appr Corr','\n\n']
tick_locations = np.arange(9)
new_labels = [ ''.join(x) for x in zip(tick_labels[0::2], tick_labels[1::2]) ]
plt.xticks(tick_locations, new_labels)
plt.ylabel('APE MD = ( abs (Corr MD - GT MD) / GT MD ) * 100 (%)')
plt.ylim([-5,50])
plt.title('Effect of LR in different tissue at different ACQ when SNR = 100')
plt.legend(loc='upper left')
plt.show()

