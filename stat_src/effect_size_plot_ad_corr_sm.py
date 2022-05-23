from numpy.random import seed
import nibabel as nib
from numpy import mean
from numpy import var
from math import sqrt
import numpy, pandas, researchpy
import numpy as np
import pandas as pd

def _cohend(true, corpt):
    cohen_d = np.zeros((112,112,54))
    for x in range(112):
        for y in range(112):
            for z in range(54):
                d1 = []
                d2 = []
                for i in range(19):
                    d1.append(true[i][x,y,z])
                    d2.append(corpt[i][x,y,z])
                n1, n2 = len(d1), len(d2)
                s1, s2 = var(d1, ddof=1), var(d2, ddof=1)
                s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
                u1, u2 = mean(d1), mean(d2)
                cohen_d[x,y,z] =  (u1 - u2) / s
    return cohen_d

def cohend(true, corpt):
    cohen_d = np.zeros((112,112,54))
    for x in range(112):
        for y in range(112):
            for z in range(54):
                d1 = []
                d2 = []
                for i in range(19):
                    d1.append(true[i][x,y,z])
                    d2.append(corpt[i][x,y,z])
                p1 = pd.Series(d1)
                p2 = pd.Series(d2)
                res = researchpy.ttest(p1, p2)
                cohen_d[x,y,z] = res[1]['results'][6]
    return cohen_d

def cohend1(true, corpt,mask):
    cohen_d = np.zeros((112,112,54))
    for x in range(112):
        for y in range(112):
            for z in range(54):
                if mask[x,y,z] == 1.0:
                    d1 = []
                    d2 = []
                    d1 = true[:,x,y,z]
                    d2 = corpt[:,x,y,z]
                    p1 = pd.Series(d1)
                    p2 = pd.Series(d2)
                    res = researchpy.ttest(p1, p2)
                    cohen_d[x,y,z] = res[1]['results'][6]
    return cohen_d

corpt_fa_dict = {}
cpre = '/nfs/masi/kanakap/projects/LR/aggregate_study/'
cpost = '/approx_corrected_ad.nii'
for i in range(1,10):
    outdir = cpre + 'OUTPUT_masivar_SNR30_' +str(i) + cpost
    corpt_fa = nib.load(outdir).get_fdata()
    corpt_fa[np.isnan(corpt_fa)] = 0
    corpt_fa_dict[i] = corpt_fa

corpt_fa_list = np.array(list(corpt_fa_dict.values()))

true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__ad.nii').get_fdata()
true_fa_dict = {}
for i in range(1,10):
    true_fa_dict[i] = true_fa
    true_fa[np.isnan(true_fa)] = 0
true_fa_list = np.array(list(true_fa_dict.values()))

mask =  nib.load('/home/local/VANDERBILT/kanakap/mask.nii').get_fdata()

aa = cohend1(true_fa_list, corpt_fa_list, mask)

img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/true__ad.nii')
coh = nib.Nifti1Image(aa, img.affine, img.header)
nib.save(coh, 'cohend_image_snr30_appr_ad.nii')
