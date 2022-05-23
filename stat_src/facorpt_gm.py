import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import operator as op
import xml.etree.cElementTree as et
import pandas as pd
import scipy.io as sio

MD_Lest = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_study/Lest_fa.nii').get_fdata()
MD_true = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()
atlas_img = nib.load('/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/slantatlas2subj.nii.gz')
atlas = atlas_img.get_fdata()
LR = sio.loadmat('../src/LRfield_posA.mat')
vL = LR['vL']

MD_diff = {}
alltrue = []
allalldiff = []
alldiff = []
L_det_roi = []
allL_det = {}
for i in range(1,208):
    for x in range(96):
        for y in range(96):
            for z in range(68):
                 if atlas[x,y,z] == i:
                     diff = (MD_Lest[x,y,z] - MD_true[x,y,z])
                     alldiff.append(diff)
                     alltrue.append(MD_true[x,y,z])
                     allalldiff.append(diff)
                     # LR field
                     L_mat = np.squeeze(vL[:,:,x,y,z]);
                     L_det_roi.append(np.linalg.det(L_mat));

    if len(alldiff) != 0:
        MD_diff[i] = alldiff
        alldiff = []
        allL_det[i] = L_det_roi
        L_det_roi=[]
print(np.nanmean(allalldiff))
# get the avg of MD diff and change the label no to that
avg_md_diff_labels = {}
for k,v in MD_diff.items():
    # v is the list of grades for student k
    avg_md_diff_labels[k] = np.nansum(v)/ float(len(v))
MD_diff_atlas = atlas.copy()
MD_diff_atlas[MD_diff_atlas == 0.0] = np.nan
for i in avg_md_diff_labels.keys():
    MD_diff_atlas[MD_diff_atlas == i] = avg_md_diff_labels[i]
nib.save(nib.Nifti1Image(MD_diff_atlas,atlas_img.affine),'/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_corpt/FAdiff_avg_gmlabels.nii.gz')

# change key to roi names
filename = '/nfs/masi/hansencb/10_29_2019_human_repositioned/3tb/posA/slant/OUTPUTS/FinalVolTxt/T1_label_volumes.txt'
with open(filename,"r") as f:
    lines=f.readlines()
    slant_roi=[]
    for x in lines:
        slant_roi.append(x.split(',')[0])
slant_roi.pop(0)

key_list = []
for i in MD_diff.keys():
    key_list.append(i)

for i,j in zip(key_list,slant_roi):
    MD_diff[j] = MD_diff[i]
    del MD_diff[i]
    allL_det[j] = allL_det[i]
    del allL_det[i]

# get the mean
Ldet_mean = {key: np.mean(allL_det[key]) for key in allL_det}
sorted_Ldet_mean = dict(sorted(Ldet_mean.items(), key=lambda item: item[1]))
sorted_MDdiff_Ldet = dict(sorted(MD_diff.items(), key=lambda kv: sorted_Ldet_mean[kv[0]]))
print(sorted_Ldet_mean.values())

# dict to df to plot in seaborn
df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in MD_diff.items() ]))
sorted_index = df.median().sort_values().index
df_sorted=df[sorted_index]

#labels, data = MD_diff.keys(), MD_diff.values()
#sns.set(rc={'figure.figsize':(11.7,50.27)})
plt.figure(num=1,figsize=(40,40))
ax = sns.boxplot(data=df_sorted)
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.ylabel('∆ FA')
#plt.xticks(range(1, len(labels) + 1), labels, rotation=90)
plt.title('Corruption of GM regions of MR scan at isocenter')

plt.figure(num=2,figsize=(40,40))
df_md_ldet = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in sorted_MDdiff_Ldet.items() ]))
#md_ldet_labels = list(df_md_ldet.columns)
ax = sns.boxplot(data=df_md_ldet)
ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
plt.ylabel('∆ FA')
#plt.yticks(range(1, len(md_ldet_labels) + 1), md_ldet_labels)
#plt.yticks(md_ldet_labels)
plt.title('Corruption of GM regions of MR scan at isocenter - sorted by LRfield')

#plt.figure(2)
#plt.scatter(alltrue,allalldiff)
#plt.ylabel('∆ FA')
#plt.xlabel('True FA')
plt.show()
