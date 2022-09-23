import os
import sys
import shutil
import tempfile
import numpy as np
import nibabel as nib
from utils import reconstruct_signal_at_voxel
from joblib import Parallel, delayed
from dipy.io import read_bvals_bvecs

# Assign inputs
n = nib.load(sys.argv[1]).get_fdata()
vol = nib.load(sys.argv[1])
og_file = sys.argv[2]
ob_file = sys.argv[3]
vec_folder = sys.argv[4]
val_folder = sys.argv[5]
out_dir = sys.argv[6]
out_prefix = sys.argv[7]

# get the shell index - to separate the shell 
og_bval, og_bvec = read_bvals_bvecs(ob_file,og_file)
ind_1000 = np.where(og_bval == 1000)
ind_2000 = np.where(og_bval == 2000)
ind_b0 = np.nonzero(og_bval == 0)
ind_b0 = np.squeeze(ind_b0)
ind_0_1000 = np.where((og_bval == 0) | (og_bval == 1000))
ind_0_2000 = np.where((og_bval == 0) | (og_bval == 2000))

# load the voxelwise bvec and bvals
bvec_vols = []
for i in sorted(os.listdir(vec_folder)):
    if i.endswith('.nii.gz'):
        bvec_vol = nib.load(vec_folder + '/' + i).get_fdata()
        bvec_vol = np.expand_dims(bvec_vol,4)
        bvec_vol = np.transpose(bvec_vol,(0,1,2,4,3))
        bvec_vols.append(bvec_vol)
bvec_stack = np.stack(bvec_vols,3)
bvec_stack = bvec_stack.squeeze()

bval_vols = []
for i in sorted(os.listdir(val_folder)):
    if i.endswith('.nii.gz'):
        bval_vol = nib.load(val_folder + '/' + i).get_fdata()
        bval_vols.append(bval_vol)
bval_stack = np.stack(bval_vols,3)

num_cores = 10
path = tempfile.mkdtemp()
xaxis = range(n.shape[0])
yaxis = range(n.shape[1])
zaxis = range(n.shape[2]) 

# for dwi with 0 1000
len1 = ind_0_1000[0]
dwi_hat_path1 = os.path.join(path,'dwi_hat1.mmap')
dwi_hat1 = np.memmap(dwi_hat_path1, dtype=float, shape=(n.shape[0],n.shape[1],n.shape[2],ind_1000[0].shape[0]), mode='w+')
data = n[:,:,:,len1]
org_bvec = og_bvec[len1,:]
org_bval = og_bval[len1]
corr_bvec = bvec_stack[:,:,:,len1,:]
corr_bval = bval_stack[:,:,:,len1]
sh_order = 6
results = Parallel(n_jobs=num_cores)(delayed(reconstruct_signal_at_voxel)(i,j,k,data,org_bvec,org_bval,corr_bvec,corr_bval,dwi_hat1,sh_order) for k in zaxis for j in yaxis for i in xaxis)

# for dwi with 0 2000
len2 = ind_0_2000[0]
dwi_hat_path2 = os.path.join(path,'dwi_hat2.mmap')
dwi_hat2 = np.memmap(dwi_hat_path2, dtype=float, shape=(n.shape[0],n.shape[1],n.shape[2],ind_2000[0].shape[0]), mode='w+')
data = n[:,:,:,len2]
org_bvec = og_bvec[len2,:]
org_bval = og_bval[len2]
corr_bvec = bvec_stack[:,:,:,len2,:]
corr_bval = bval_stack[:,:,:,len2]
sh_order = 8
results = Parallel(n_jobs=num_cores)(delayed(reconstruct_signal_at_voxel)(i,j,k,data,org_bvec,org_bval,corr_bvec,corr_bval,dwi_hat2,sh_order) for k in zaxis for j in yaxis for i in xaxis)

dwmri_corrected = np.zeros((n.shape[0],n.shape[1],n.shape[2],n.shape[3]))

dwmri_corrected[:,:,:,ind_b0] = n[:,:,:,ind_b0] 
dwmri_corrected[:,:,:,ind_1000[0]] = dwi_hat1
dwmri_corrected[:,:,:,ind_2000[0]] = dwi_hat2
dwmri_corrected = np.nan_to_num(dwmri_corrected)

out_name = out_prefix + '_Lcorrected_sig.nii.gz'
output_file = os.path.join(out_dir, out_name )

nib.save(nib.Nifti1Image(dwmri_corrected.astype(np.float32),vol.affine),output_file )


# Remove tmp file
try:
    shutil.rmtree(path)
except:
    print("Couldn't delete folder")


