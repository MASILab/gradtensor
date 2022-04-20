import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import nmmn.plots
parula=nmmn.plots.parulacmap() # for MATLAB's cmap
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__fa.nii').get_fdata()

corpt_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/uncorrected_fa.nii').get_fdata()
sm_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/approx_corrected_fa.nii').get_fdata()
bx_corr_fa = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/emp_corrected_fa.nii').get_fdata()

corpt_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/uncorrected_fa.nii').get_fdata()
sm_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/approx_corrected_fa.nii').get_fdata()
bx_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/emp_corrected_fa.nii').get_fdata()

corpt_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/uncorrected_fa.nii').get_fdata()
sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/approx_corrected_fa.nii').get_fdata()
bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1/emp_corrected_fa.nii').get_fdata()

corpt_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/uncorrected_fa.nii').get_fdata()
sm_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/approx_corrected_fa.nii').get_fdata()
bx_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1/emp_corrected_fa.nii').get_fdata()

seg =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1/seg_reg_skull_stripped.nii.gz').get_fdata()

# wm
def plot_one_snr(corpt_fa,true_fa,snr):
    err_fa = corpt_fa - true_fa
    pe_fa = 100 * err_fa / true_fa
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa_csf = np.where(seg == 0, ape_fa , math.nan) # for csf 
    ape_fa_wm = np.where(seg == 1, ape_fa , math.nan) 
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
    cdf = pd.concat([df1, df2, df3])
    cdf = cdf.assign(SNR = snr)
    return cdf

cdf = plot_one_snr(corpt_fa,true_fa,'inf')
cdf_ec = plot_one_snr(bx_corr_fa,true_fa,'inf ec')

cdf100 = plot_one_snr(corpt_fa100,true_fa,'100')
cdf30 = plot_one_snr(corpt_fa30,true_fa,'30')
cdf10 = plot_one_snr(corpt_fa10,true_fa,'10')

all_cdf = pd.concat([cdf,cdf_ec, cdf100, cdf30, cdf10])
all_cdf = all_cdf.rename(columns={0:'Abs percent error (%)'})
plt.figure()
sns.set(font_scale = 2)
ax = sns.boxplot(data=all_cdf,hue='Tissue',x = 'SNR',y='Abs percent error (%)')
ax.set_xticklabels(['','inf','100','30','10'])
plt.ylim([-0.5,500])
plt.title('Effect of LR in different tissue at different SNR')
plt.show()

"""
err_fa = bx_corr_fa - true_fa
pe_fa = 100 * err_fa / true_fa
ape_fa = np.abs(pe_fa)
ape_fa[np.isnan(ape_fa)] = 0
ape_fa_csf = np.where(seg == 0, ape_fa , math.nan) # for csf
ape_fa_wm = np.where(seg == 1, ape_fa , math.nan)
ape_fa_gm = np.where(seg == 2, ape_fa , math.nan)
ape_fa_csf = ape_fa_csf.reshape(-1)
ape_fa_csf = [x for x in ape_fa_csf if math.isnan(x) == False]
ape_fa_gm = ape_fa_gm.reshape(-1)
ape_fa_gm = [x for x in ape_fa_gm if math.isnan(x) == False]
ape_fa_wm = ape_fa_wm.reshape(-1)
ape_fa_wm = [x for x in ape_fa_wm if math.isnan(x) == False]
df1 = pd.DataFrame(ape_fa_gm).assign(Tissue='GM')
df2 = pd.DataFrame(ape_fa_wm).assign(Tissue='WM')
df3 = pd.DataFrame(ape_fa_csf).assign(Tissue='CSF')
cdf_ec = pd.concat([df1, df2, df3])
cdf_ec = cdf_ec.assign(SNR = 'inf EC')

err_fa100 = corpt_fa100 - true_fa
pe_fa100 = 100 * err_fa100 / true_fa
ape_fa100 = np.abs(pe_fa100)
ape_fa_csf100 = np.where(seg == 0, ape_fa100 ,math.nan) # for csf
ape_fa_wm100 = np.where(seg == 1, ape_fa100 ,math.nan)
ape_fa_gm100 = np.where(seg == 2, ape_fa100 ,math.nan)
ape_fa_csf100 = ape_fa_csf100.reshape(-1)
ape_fa_gm100 = ape_fa_gm100.reshape(-1)
ape_fa_wm100 = ape_fa_wm100.reshape(-1)
ape_fa_csf100 = [x for x in ape_fa_csf100 if math.isnan(x) == False]
ape_fa_gm100 = [x for x in ape_fa_gm100 if math.isnan(x) == False]
ape_fa_wm100 = [x for x in ape_fa_wm100 if math.isnan(x) == False]
df1 = pd.DataFrame(ape_fa_gm100).assign(Tissue='GM')
df2 = pd.DataFrame(ape_fa_wm100).assign(Tissue='WM')
df3 = pd.DataFrame(ape_fa_csf100).assign(Tissue='CSF')
cdf100 = pd.concat([df1, df2, df3])
cdf100 = cdf100.assign(SNR = '100')

err_fa30 = corpt_fa30 - true_fa
pe_fa30 = 100 * err_fa30 / true_fa
ape_fa30 = np.abs(pe_fa30)
ape_fa_csf30 = np.where(seg == 0, ape_fa30 ,math.nan) # for csf
ape_fa_wm30 = np.where(seg == 1, ape_fa30 ,math.nan)
ape_fa_gm30 = np.where(seg == 2, ape_fa30 ,math.nan)
ape_fa_csf30 = ape_fa_csf30.reshape(-1)
ape_fa_gm30 = ape_fa_gm30.reshape(-1)
ape_fa_wm30 = ape_fa_wm30.reshape(-1)
ape_fa_csf30 = [x for x in ape_fa_csf30 if math.isnan(x) == False]
ape_fa_gm30 = [x for x in ape_fa_gm30 if math.isnan(x) == False]
ape_fa_wm30 = [x for x in ape_fa_wm30 if math.isnan(x) == False]
df1 = pd.DataFrame(ape_fa_gm30).assign(Tissue='GM')
df2 = pd.DataFrame(ape_fa_wm30).assign(Tissue='WM')
df3 = pd.DataFrame(ape_fa_csf30).assign(Tissue='CSF')
cdf30 = pd.concat([df1, df2, df3])
cdf30 = cdf30.assign(SNR = '30')

err_fa10 = corpt_fa10 - true_fa
pe_fa10 = 100 * err_fa10 / true_fa
ape_fa10 = np.abs(pe_fa10)
ape_fa_csf10 = np.where(seg == 0, ape_fa10 ,math.nan) # for csf
ape_fa_wm10 = np.where(seg == 1, ape_fa10 ,math.nan)
ape_fa_gm10 = np.where(seg == 2, ape_fa10 ,math.nan)
ape_fa_csf10 = ape_fa_csf10.reshape(-1)
ape_fa_gm10 = ape_fa_gm10.reshape(-1)
ape_fa_wm10 = ape_fa_wm10.reshape(-1)
ape_fa_csf10 = [x for x in ape_fa_csf10 if math.isnan(x) == False]
ape_fa_gm10 = [x for x in ape_fa_gm10 if math.isnan(x) == False]
ape_fa_wm10 = [x for x in ape_fa_wm10 if math.isnan(x) == False]
df1 = pd.DataFrame(ape_fa_gm10).assign(Tissue='GM')
df2 = pd.DataFrame(ape_fa_wm10).assign(Tissue='WM')
df3 = pd.DataFrame(ape_fa_csf10).assign(Tissue='CSF')
cdf10 = pd.concat([df1, df2, df3])
cdf10 = cdf10.assign(SNR = '10')

all_cdf = pd.concat([cdf,cdf_ec, cdf100, cdf30, cdf10])
all_cdf = all_cdf.rename(columns={0:'Abs percent error (%)'})
plt.figure()
sns.set(font_scale = 2)
ax = sns.boxplot(data=all_cdf,hue='Tissue',x = 'SNR',y='Abs percent error (%)')
ax.set_xticklabels(['','inf','100','30','10'])
plt.ylim([-0.5,500])
plt.title('Effect of LR in different tissue at different SNR')
plt.show()
"""
