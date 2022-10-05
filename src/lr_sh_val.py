import os
import sys
import math
import shutil
import tempfile
import numpy as np
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as plt
from dipy.io import read_bvals_bvecs
from joblib import Parallel, delayed
from scilpy.reconst.sh import compute_rish
from utils import val_pk, val_emp, angularCorrCoeff
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Assign inputs
n = nib.load(sys.argv[1]).get_fdata()
vol = nib.load(sys.argv[1])
og_file = sys.argv[2]
ob_file = sys.argv[3]
vec_folder = sys.argv[4]
val_folder = sys.argv[5]
emp_sig = nib.load(sys.argv[6]).get_fdata()
mask = nib.load(sys.argv[7]).get_fdata()

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
sh_acc1000 = sh_acc1000[mask==1] 
sh_acc2000 = sh_acc2000[mask==1]

print('ACC SH 1000 ',np.nanmean(np.nanmean(sh_acc1000)))
print('ACC SH 2000 ',np.nanmean(np.nanmean(sh_acc2000)))

# RISH features

emp_rish1000, final_orders = compute_rish(emp_sh1000,mask=mask)
pk_rish1000, final_orders = compute_rish(pk_sh1000,mask=mask)

emp_rish2000, final_orders = compute_rish(emp_sh2000,mask=mask)
pk_rish2000, final_orders = compute_rish(pk_sh2000,mask=mask)

rmse1000 = []
rmse2000 = []
corr1000 = []
corr2000 = []
pvalue1000 = []
pvalue2000 = []

for i in range(4):
    a = emp_rish1000[:,:,:,i]
    b = pk_rish1000[:,:,:,i]
    MSE = np.square(np.subtract(a,b)).mean() 
    RMSE = math.sqrt(MSE)
    rmse1000.append(RMSE)
    ra = a.flatten()
    rb = b.flatten()
    corr, _ = stats.pearsonr(ra,rb)
    corr1000.append(corr)
    p = stats.wilcoxon(ra,rb)
    pvalue1000.append(p.pvalue)

for i in range(5):
    a = emp_rish2000[:,:,:,i]
    b = pk_rish2000[:,:,:,i]
    MSE = np.square(np.subtract(a,b)).mean() 
    RMSE = math.sqrt(MSE)
    rmse2000.append(RMSE)
    ra = a.flatten()
    rb = b.flatten()
    corr, _ = stats.pearsonr(ra,rb)
    corr2000.append(corr)
    p = stats.wilcoxon(ra,rb)
    pvalue2000.append(p.pvalue)


print('ACC RISH 0 1000 ',rmse1000)
print('ACC RISH 0 2000 ',rmse2000)

print('P CORR RISH 0 1000 ',corr1000)
print('P CORR RISH 0 2000 ',corr2000)

print('P VALUE RISH 0 1000 ',pvalue1000)
print('P VALUE RISH 0 2000 ',pvalue2000)

def plot_rish(features, ax, row, col):
    slice_idx = 60
    slice = features
    slice = slice[slice_idx, :,:]
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    ax[row,col].axis('off')
    cmap = plt.get_cmap('viridis')
    m = 0
    M = 5
    im = ax[row,col].imshow(np.abs(slice), vmin=m, vmax=M, cmap=cmap)
    divider = make_axes_locatable(ax[row,col])
    ax = divider.append_axes("right", size="5%", pad=0.05)
    a = plt.colorbar(im, cax=ax)

plt.rcParams.update({'font.size':12})
fig, ax = plt.subplots(2,3,figsize=(10,3))
plot_rish(emp_rish1000[:,:,:,0] , ax, 0,0)
plot_rish(pk_rish1000[:,:,:,0] , ax, 0,1)
plot_rish((emp_rish1000[:,:,:,0] - pk_rish1000[:,:,:,0]), ax, 0,2)

plot_rish(emp_rish2000[:,:,:,0] , ax,1,0)
plot_rish(pk_rish2000[:,:,:,0] , ax,1,1)
plot_rish((emp_rish2000[:,:,:,0] - pk_rish2000[:,:,:,0]),ax, 1,2)
plt.savefig('rish_difference.png')

try:
    shutil.rmtree(path)
except:
    print("Couldn't delete folder")
