import os
import sys
from glob import glob
import nibabel as nib
import numpy as np


#a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_est/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()
#b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/sub-cIVs001_ses-s1Bx2__'+sys.argv[1]+'.nii.gz').get_fdata()

#use -> a = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs052/ses-s1Bx1/tracto_op_lr_corr_1/sub-cIVs052_ses-s1Bx1__'+sys.argv[1]+'.nii.gz').get_fdata()
#use -> b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs052/ses-s1Bx1/tracto_op_corpt_1/sub-cIVs052_ses-s1Bx1__'+sys.argv[1]+'.nii.gz').get_fdata()

#a =  nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/true_'+sys.argv[1]+'.nii').get_fdata()
#b = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_1/corpt_'+sys.argv[1]+'.nii').get_fdata()

def angular_error(PEa, PEb, halfPi=True):
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
            np.square(PEa[..., 1] - PEb[..., 1]) + \
            np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)

    ang = 2 * np.real(np.arcsin(chord/2))

    if halfPi:
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]
    return np.degrees(ang)

#mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/bundles/*_mask.nii.gz'))
#mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1_Lest/tractsegout/TOM_trackings/*_mask.nii.gz'))
#mask_files = glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs052/ses-s1Bx1/tracto_op_lr_corr_1/bundles/*_mask.nii.gz'))

# get fa files for only session in each subject
p = []
fa_corr_files = []
fa_files = []
mask_folder = []
for path in glob(os.path.join('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/*/')):
        p.append(path)
p.sort()
for i in range(len(p)):
    count = 0
    for j in os.listdir(p[i]):
        count = count+1
        if count != 2:
            f = p[i] + j + '/tracto_op_corpt_1/'+p[i].split('/')[7] + '_' + j + '__evecs_v1.nii.gz'
            fa_files.append(f)
for i in range(len(p)):
    count = 0
    for j in os.listdir(p[i]):
        count = count+1
        if count != 2:
            f = p[i] + j + '/tracto_op_lr_corr_1/'+p[i].split('/')[7] + '_' + j + '__evecs_v1.nii.gz'
            fa_corr_files.append(f)
for i in range(len(p)):
    count = 0
    for j in os.listdir(p[i]):
        count = count+1
        if count != 2:
            f = p[i] + j 
            mask_folder.append(f)
mask_folder.sort()
fa_corr_files.sort()
fa_files.sort()

#ang_error= angular_error(b,a,halfPi=True)
values = []
bundle_ape = {}
ae =[]
allvalue = {}
for j in range(len(fa_files)):
    print(j)
    a = nib.load(fa_corr_files[j]).get_fdata()
    b = nib.load(fa_files[j]).get_fdata()
    mask_files = glob(os.path.join(mask_folder[j]+'/tracto_op_lr_corr_1/bundles/*_mask.nii.gz'))
    ang_error= angular_error(b,a,halfPi=True)
    for i in range(len(mask_files)):
        mask = nib.load(mask_files[i]).get_fdata()
        for x in range(a.shape[0]):
            for y in range(a.shape[1]):
                for z in range(a.shape[2]):
                    if mask[x,y,z] == 1:
                        err_fa = ang_error[x,y,z]
                        ae.append(err_fa)
        values.append(np.nanmedian(ae)) # the values for avg ape in each bundle
        bundle_name = mask_files[i].split('/')[11].split('__')[1].split('_clean')[0]
        bundle_ae[bundle_name] = np.nanmedian(ae)
        ae = []
    allvalue[j] = bundle_ae
    bundle_ae = {}
df = pd.DataFrame.from_dict(allvalue)
mean = df.mean(axis=1)
print(values)
print(min(values))
print(max(values))
print(mask_files[values.index(max(values))]) # maximum mask bundle
