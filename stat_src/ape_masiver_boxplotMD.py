import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import nmmn.plots
#parula=nmmn.plots.parulacmap() # for MATLAB's cmap
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__md.nii').get_fdata()

corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/uncorrected_md.nii').get_fdata()
sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/approx_corrected_md.nii').get_fdata()
bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/emp_corrected_md.nii').get_fdata()

corpt_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_md.nii').get_fdata()
sm_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/emp_corrected_md.nii').get_fdata()

corpt_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/uncorrected_md.nii').get_fdata()
sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/emp_corrected_md.nii').get_fdata()

corpt_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/uncorrected_md.nii').get_fdata()
sm_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/approx_corrected_md.nii').get_fdata()
bx_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/emp_corrected_md.nii').get_fdata()

seg =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/seg_reg_skull_stripped.nii.gz').get_fdata()

# wm
def plot_one_snr(corpt_fa,true_fa,snr):
    err_fa = corpt_fa - true_fa
    ape_fa =  (abs(err_fa) / true_fa) * 100
    #ape_fa = np.abs(pe_fa)
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
    cdf = cdf.assign(SNR = snr)
    median = cdf.groupby(['SNR', 'Tissue'])[0].median().values
    idx = [0,1,2]
    median_sort = median[idx]
    return cdf, median

cdf,median = plot_one_snr(corpt_fa,true_fa,'inf')
cdf_ec,median_ec = plot_one_snr(bx_corr_fa,true_fa,'inf ec')
cdf_ac,median_ac = plot_one_snr(sm_corr_fa,true_fa,'inf ac')

cdf100,median100 = plot_one_snr(corpt_fa100,true_fa,'100')
cdf_ec100,median_ec100 = plot_one_snr(bx_corr_fa100,true_fa,'100 ec')
cdf_ac100,median_ac100 = plot_one_snr(sm_corr_fa100,true_fa,'100 ac')

cdf30,median30 = plot_one_snr(corpt_fa30,true_fa,'30')
cdf_ec30,median_ec30 = plot_one_snr(bx_corr_fa30,true_fa,'30 ec')
cdf_ac30,median_ac30 = plot_one_snr(sm_corr_fa30,true_fa,'30 ac')

cdf10,median10 = plot_one_snr(corpt_fa10,true_fa,'10')
cdf_ec10,median_ec10 = plot_one_snr(bx_corr_fa10,true_fa,'10 ec')
cdf_ac10,median_ac10 = plot_one_snr(sm_corr_fa10,true_fa,'10 ac')

all_cdf = pd.concat([cdf,cdf_ec,cdf_ac,cdf100,cdf_ec100,cdf_ac100, cdf30, cdf_ec30,cdf_ac30,cdf10,cdf_ec10,cdf_ac10])
all_cdf = all_cdf.rename(columns={0:'Abs percent error (%)'})
plt.figure()
sns.set(font_scale = 1.3)
flierprops = dict(markerfacecolor='0.75', markersize=2,
              linestyle='none')
testPlot = sns.boxplot(data=all_cdf,hue='Tissue',x = 'SNR',y='Abs percent error (%)',palette="Set2",flierprops=flierprops)

#m1 = all_cdf.groupby(['SNR', 'Tissue'])['Abs percent error (%)'].median().values
m1 = np.concatenate([median,median_ec,median_ac,median100,median_ec100,median_ac100,median30,median_ec30,median_ac30,median10,median_ec10,median_ac10])
mL1 = [str(np.round(s, 2)) for s in m1]
ind = 0
for tick in range(len(testPlot.get_xticklabels())):
    testPlot.text(tick+.4, m1[ind+2]+1, mL1[ind+2],  horizontalalignment='center',  color='k', weight='semibold', fontsize=10)
    testPlot.text(tick+.1, m1[ind+1]+1, mL1[ind+1],  horizontalalignment='center',  color='k', weight='semibold', fontsize=10)
    testPlot.text(tick-.3, m1[ind]+1, mL1[ind], horizontalalignment='center', color='k', weight='semibold', fontsize=10)
    ind += 3 

#testPlot.set_xticklabels(['','inf','','','100','','','30','','','10',''])
#testPlot.set_xticklabels(['','inf','','','100','','','30','','','10',''])
#tick_labels=['','\n\n ','Corruption','\n\nFri',' ','\n\nSat',' ','\n\nSun','Emp Corr','\n\nMon',' ','\n\nTue',' ','\n\nWed','Appr Corr',
 #            '\n\nThu',' ','\n\nFri',' ','\n\nSat','Corrpution','\n\nSun',' ','\n\nMon',' ','\n\nTue','Emp Corr','\n\nWed',' ',
 #            '\n\nThu',' ','\n\nFri',' ','\n\nSat','Appr Corr','\n\nSun',' ','\n\nMon',' ','\n\nTue','Corrpution','\n\nWed',' ',
 #            '\n\nThu',' ','\n\nFri','Emp Corr','\n\nSat',' ','\n\nSun',' ','\n\nMon','Appr Corr','\n\nTue',' ','\n\nWed']
# tick_labels=['Corruption','\n\n ','Emp Corr','\n\ninf','Appr Corr','\n\n','Corruption','Emp Corr','Appr Corr','Corruption','Emp Corr','Appr Corr','Corruption','Emp Corr','Appr Corr']
tick_labels=['Corruption','\n\n ','Emp Corr','\n\ninf','Appr Corr','\n\n','Corruption','\n\n ','Emp Corr','\n\n100','Appr Corr','\n\n','Corruption','\n\n ','Emp Corr','\n\n30','Appr Corr','\n\n','Corruption','\n\n ','Emp Corr','\n\n10','Appr Corr','\n\n']

tick_locations = np.arange(12)

new_labels = [ ''.join(x) for x in zip(tick_labels[0::2], tick_labels[1::2]) ]
#plt.figure(figsize=(15, 8))

plt.xticks(tick_locations, new_labels)

plt.ylim([-10,75])
plt.title('Effect of LR in different tissue at different SNR')
plt.ylabel('APE MD = ( abs (Corr MD - GT MD) / GT MD ) * 100 (%)')
plt.legend(loc='upper left')
plt.show()

