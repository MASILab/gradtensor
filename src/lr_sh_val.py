import os
import sys
import math
import shutil
import tempfile
import numpy as np
import nibabel as nib
from dipy.io import read_bvals_bvecs
from joblib import Parallel, delayed
from scilpy.reconst.sh import compute_rish
from utils import val_pk, val_emp, angularCorrCoeff


# Assign inputs
n = nib.load(sys.argv[1]).get_fdata()
vol = nib.load(sys.argv[1])
og_file = sys.argv[2]
ob_file = sys.argv[3]
vec_folder = sys.argv[4]
val_folder = sys.argv[5]
emp_sig = nib.load(sys.argv[6]).get_fdata()

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

# Get sh for pk emp sig
# for 1000
len1 = ind_0_1000[0]
sh_order = 6
pk_sh1000 = val_pk(emp_sig[:,:,:,len1],og_bval[len1],og_bvec[len1,:],sh_order)
# for 2000
len2 = ind_0_2000[0]
sh_order = 8
pk_sh2000 = val_pk(emp_sig[:,:,:,len2],og_bval[len2],og_bvec[len2,:],sh_order)

# Get sh for emp vox bvec bval
num_cores = 10
path = tempfile.mkdtemp()
xaxis = range(n.shape[0])
yaxis = range(n.shape[1])
zaxis = range(n.shape[2]) 

# for dwi with 0 1000
len1 = ind_0_1000[0]
dwi_hat_path1 = os.path.join(path,'emp_sh1000.mmap')
emp_sh1000 = np.memmap(dwi_hat_path1, dtype=float, shape=(n.shape[0],n.shape[1],n.shape[2],28), mode='w+')
data = n[:,:,:,len1]
org_bval = og_bval[len1]
corr_bvec = bvec_stack[:,:,:,len1,:]
corr_bval = bval_stack[:,:,:,len1]
sh_order = 6
results = Parallel(n_jobs=num_cores)(delayed(val_emp)(i,j,k,data,org_bval,corr_bvec,corr_bval,emp_sh1000,sh_order) for k in zaxis for j in yaxis for i in xaxis)

# for dwi with 0 2000
len2 = ind_0_2000[0]
dwi_hat_path2 = os.path.join(path,'emp_sh2000.mmap')
emp_sh2000 = np.memmap(dwi_hat_path2, dtype=float, shape=(n.shape[0],n.shape[1],n.shape[2],45), mode='w+')
data = n[:,:,:,len2]
org_bval = og_bval[len2]
corr_bvec = bvec_stack[:,:,:,len2,:]
corr_bval = bval_stack[:,:,:,len2]
sh_order = 8
results = Parallel(n_jobs=num_cores)(delayed(val_emp)(i,j,k,data,org_bval,corr_bvec,corr_bval,emp_sh2000,sh_order) for k in zaxis for j in yaxis for i in xaxis)

sh_acc1000 = angularCorrCoeff(emp_sh1000,pk_sh1000)
sh_acc2000 = angularCorrCoeff(emp_sh2000,pk_sh2000)

print('ACC SH 1000 ',np.nanmean(np.nanmean(sh_acc1000)))
print('ACC SH 2000 ',np.nanmean(np.nanmean(sh_acc2000)))

# RISH features

emp_rish1000, final_orders = compute_rish(emp_sh1000)
pk_rish1000, final_orders = compute_rish(pk_sh1000)

emp_rish2000, final_orders = compute_rish(emp_sh2000)
pk_rish2000, final_orders = compute_rish(pk_sh2000)

rmse1000 = []
rmse2000 = []

for i in range(4):
    MSE = np.square(np.subtract(emp_rish2000[:,:,:,i],pk_rish2000[:,:,:,i])).mean() 
    RMSE = math.sqrt(MSE)
    rmse1000.append(np.nanmean(np.nanmean(RMSE)))

for i in range(5):
    MSE = np.square(np.subtract(emp_rish2000[:,:,:,i],pk_rish2000[:,:,:,i])).mean() 
    RMSE = math.sqrt(MSE)
    rmse2000.append(RMSE)

print('ACC RISH 0 1000 ',rmse1000)
print('ACC RISH 0 2000 ',rmse2000)


try:
    shutil.rmtree(path)
except:
    print("Couldn't delete folder")