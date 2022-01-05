import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import operator as op
import xml.etree.cElementTree as et
import pandas as pd
import scipy.io as sio

MD_Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx/Lest_corr_bx_primary_eigvec.nii').get_fdata()
MD_true = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_primary_eigvec.nii').get_fdata()
atlas_img = nib.load('/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/reg/FAatlas2subj.nii.gz')
atlas = atlas_img.get_fdata()
LR = sio.loadmat('../src/LRfield_posA.mat')
vL = LR['vL']

def angular_error(PEa, PEb, halfPi=True):
    # PEa = preprocessing.normalize(PEa, norm='l1', axis=1)
    # PEb = preprocessing.normalize(PEb, norm='l1', axis=1)
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
            np.square(PEa[..., 1] - PEb[..., 1]) + \
            np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)

    ang = 2 * np.real(np.arcsin(chord/2))

    if halfPi:
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]

    return np.degrees(ang)

ang_error = angular_error(MD_Lest,MD_true,halfPi=True)

MD_diff = {}
alltrue = []
allalldiff = []
alldiff = []
L_det_roi = []
allL_det = {}
for i in range(1,51):
    for x in range(96):
        for y in range(96):
            for z in range(68):
                 if atlas[x,y,z] == i:
                     diff = ang_error[x,y,z]
                     alldiff.append(diff)
                     alltrue.append(MD_true[x,y,z])
                     allalldiff.append(diff)
                     # LR field
                     L_mat = np.squeeze(vL[:,:,x,y,z]);
                     L_det_roi.append(np.linalg.det(L_mat));

    #if len(alldiff) != 0:
    MD_diff[i] = alldiff
    alldiff = []
    allL_det[i] = L_det_roi
    L_det_roi=[]

# get the avg of MD diff and change the label no to that
avg_md_diff_labels = {}
for k,v in MD_diff.items():
    # v is the list of grades for student k
    avg_md_diff_labels[k] = sum(v)/ float(len(v))
MD_diff_atlas = atlas.copy()
MD_diff_atlas[MD_diff_atlas == 0.0] = np.nan
for i in range(1,51):
    MD_diff_atlas[MD_diff_atlas == i] = avg_md_diff_labels[i]
nib.save(nib.Nifti1Image(MD_diff_atlas,atlas_img.affine),'/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx/pevdiff_avg_wmlabels.nii.gz')

# change key to roi names
# change key to roi names
tree=et.parse('/home/local/VANDERBILT/kanakap/fsl_605/data/atlases/JHU-labels.xml')
root=tree.getroot()
roi = []
for i in root.iter('label'):
    roi.append(i.text)
for i in range(1,51):
    MD_diff[roi[i]] = MD_diff.pop(i)
    allL_det[roi[i]] = allL_det.pop(i)

# get the mean for ldet, sort ldet, map those keys to MD diff keys
Ldet_mean = {key: np.mean(allL_det[key]) for key in allL_det}
sorted_Ldet_mean = dict(sorted(Ldet_mean.items(), key=lambda item: item[1]))
sorted_MDdiff_Ldet = dict(sorted(MD_diff.items(), key=lambda kv: sorted_Ldet_mean[kv[0]]))
print(sorted_Ldet_mean.values())

# dict to df to plot in seaborn
df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in MD_diff.items() ]))
sorted_index = df.median().sort_values().index
df_sorted=df[sorted_index]

med_labels = list(df_sorted.columns)
plt.figure(num=1,figsize=(40,40))
sns.boxplot(data=df_sorted, orient="h")
plt.xlabel('∆ PEV degrees')
#plt.yticks(range(1, len(med_labels) + 1), med_labels)
#plt.yticks(med_labels)
plt.title('Full Correction of WM regions of MR scan at isocenter')
plt.xlim([-10, 80])
plt.subplots_adjust(left=0.445, bottom=0.065, right=0.985, top=0.950, wspace=0, hspace=0)

plt.figure(num=2,figsize=(40,40))
df_md_ldet = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in sorted_MDdiff_Ldet.items() ]))
md_ldet_labels = list(df_md_ldet.columns)
sns.boxplot(data=df_md_ldet, orient="h")
plt.xlabel('∆ PEV degrees')
#plt.yticks(range(1, len(md_ldet_labels) + 1), md_ldet_labels)
#plt.yticks(md_ldet_labels)
plt.xlim([-10, 80])
plt.title('Full Correction of WM regions of MR scan at isocenter - sorted by LRfield')
plt.subplots_adjust(left=0.445, bottom=0.065, right=0.985, top=0.950, wspace=0, hspace=0)
plt.show()



