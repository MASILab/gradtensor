from numpy import mean
from numpy import var
from math import sqrt
import numpy, pandas, researchpy
import numpy as np
import pandas as pd
import sys
import nibabel as nib 

def cohend1(true, corpt,mask,img_shape):
    cohen_d = np.zeros((img_shape[0],img_shape[1],img_shape[2]))
    for x in range(img_shape[0]):
        for y in range(img_shape[1]):
            for z in range(img_shape[2]):
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

def get_list(cpost):
    data_dict = {}
    cpre = '/nfs/masi/kanakap/projects/LR/masivar_output/'
    for i in range(1,10):
        outdir = cpre + 'SNR100_d32_' +str(i) + cpost
        data = nib.load(outdir).get_fdata()
        data[np.isnan(data)] = 0
        data_dict[i] = data
    data_list = np.array(list(data_dict.values()))
    return data_list

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+sys.argv[1]+'.nii').get_fdata()
img_shape = true_fa.shape
true_fa_dict = {}
for i in range(1,10):
    true_fa_dict[i] = true_fa
    true_fa[np.isnan(true_fa)] = 0
true_fa_list = np.array(list(true_fa_dict.values()))

mask = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/mask.nii').get_fdata()

corpt_fa_list = get_list('/uncorrected_'+sys.argv[1]+'.nii')
noise_fa_list = get_list('/uncorrected_Nest_'+sys.argv[1]+'.nii')
emp_corr_fa_list = get_list('/emp/emp_corrected_'+sys.argv[1]+'.nii')

uncorrected = cohend1(true_fa_list, corpt_fa_list, mask,img_shape)
corrected = cohend1(emp_corr_fa_list, noise_fa_list, mask,img_shape)

img = nib.load('/nfs/masi/kanakap/projects/LR/masivar_input/1/true_'+sys.argv[1]+'.nii')
uncorrected_nii = nib.Nifti1Image(uncorrected, img.affine, img.header)
nib.save(uncorrected_nii, '/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_cohen/uncorrected_snr30_'+sys.argv[1]+'.nii')
corrected_nii = nib.Nifti1Image(corrected, img.affine, img.header)
nib.save(corrected_nii, '/nfs/masi/kanakap/projects/LR/masivar_output/SNR30_d32_cohen/corrected_snr30_'+sys.argv[1]+'.nii')

