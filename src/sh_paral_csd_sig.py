import os
import scipy
import tempfile
import numpy as np
import nibabel as nib
import multiprocessing
from dipy.data import get_sphere
from dipy.io import read_bvals_bvecs
from joblib import Parallel, delayed
from dipy.core.gradients import gradient_table
from scilpy.reconst.multi_processes import fit_from_model, convert_sh_basis
from scilpy.reconst.raw_signal import compute_sh_coefficients
from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel

n = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/2000dwi_fodf.nii.gz').get_fdata()
vol = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/2000dwi_fodf.nii.gz')
mask = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/mask.nii.gz').get_fdata()
vec_folder = '/home/local/VANDERBILT/kanakap/try_fx_emp/bvec'
obvecval_folder = '/home/local/VANDERBILT/kanakap/try_fx_emp/bval'
og_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/2000bvec_fodf.bvec'
ob_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/2000bval_fodf.bval'
obval, obvec = read_bvals_bvecs(ob_file,og_file)

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

reg_sphere = get_sphere('symmetric362')
#odf = np.zeros((n.shape[0],n.shape[1],n.shape[2],45))
path = tempfile.mkdtemp()
odfpath = os.path.join(path,'odf.mmap')
odf = np.memmap(odfpath, dtype=float, shape=(n.shape[0],n.shape[1],n.shape[2],45), mode='w+')

import time
start = time.process_time()
xaxis = range(112)
yaxis = range(112)
zaxis = range(54)

def appendSpherical_np(xyz):
    ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2) # r
    ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down # theta
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0]) # pi
    return ptsnew

sphcords = appendSpherical_np(obvec)
theta = sphcords[:,4]
phi = sphcords[:,5]

def processInput(i,j,k,odf):
    if mask[i,j,k] == 1:
        print(i,j,k)
        vec = bvec_stack[i,j,k,:,:]
        val = bval_stack[i,j,k,:]
        dwi = n[i][j][k]
        dwi = np.expand_dims(dwi,(1,2,3))
        dwi = np.transpose(dwi,(1,2,3,0))
        gtab = gradient_table(val, vec)
        shm_coeff = compute_sh_coefficients(dwi, gtab, 8, 'descoteaux07')
        csd_model = ConstrainedSphericalDeconvModel(gtab, ((1.509867245705930885e-03,3.422135009646922127e-04,3.422135009646922127e-04), 1.058279199218750000e+04),reg_sphere=reg_sphere,sh_order=8)
        am = shm_coeff
        m = 8 
        n = 8
        ym = scipy.special.sph_harm(m, n, theta, phi, out=None)
        if shm_coeff.ndim == 4:
            shm_coeff = convert_sh_basis(shm_coeff, reg_sphere)
            print('inside if')
            odf[i,j,k,:] = shm_coeff

n = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/dwi_fodf.nii.gz').get_fdata()
og_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/bvec_fodf.bvec'
ob_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/bval_fodf.bval'
val, vec = read_bvals_bvecs(ob_file,og_file)
sphcords = appendSpherical_np(vec[7:13,:])
theta = sphcords[:,4]
phi = sphcords[:,5]
# for one voxel
dwi = n[i][j][k][7:13]
dwi = np.expand_dims(dwi,(1,2,3))
dwi = np.transpose(dwi,(1,2,3,0))
gtab = gradient_table(val[7:13], vec[7:13,:])
shm_coeff = compute_sh_coefficients(dwi, gtab, 2, 'descoteaux07')
order = [0,2]#,4,6,8]
for i in range(6):
     signal= []
     for bvec in range(6):
         temp = 0
         for l in order:
             for m in range(-l,l):
                 ym = scipy.special.sph_harm(m, l, theta[bvec], phi[bvec], out=None)
                 if m < 0:
                     Y = np.sqrt(2) * (-1)**m * ym.imag
                 elif m > 0:
                     Y = np.sqrt(2) * (-1)**m * ym.real
                 temp+= shm_coeff[0][0][0][i] * Y
         signal.append(temp)



     signal= []
     for bvec in range(6):
         temp = 0
         for l in order:
             for m in range(-l,l):
                 ym = scipy.special.sph_harm(m, l, theta[bvec], phi[bvec], out=None)
                 if m < 0:
                     Y = np.sqrt(2) * (-1)**m * ym.imag
                 elif m > 0:
                     Y = np.sqrt(2) * (-1)**m * ym.real
                 for i in range(6):
                     temp+= shm_coeff[0][0][0][i] * Y
         signal.append(temp)

num_cores = 10
#results = Parallel(n_jobs=num_cores,require='sharedmem',mmap_mode='w+',max_nbytes='500M')(delayed(processInput)(i,j,k,odf) for k in zaxis for j in yaxis for i in xaxis)
results = Parallel(n_jobs=num_cores)(delayed(processInput)(i,j,k,odf) for k in zaxis for j in yaxis for i in xaxis)
print(time.process_time() - start)

nib.save(nib.Nifti1Image(odf.astype(np.float32),vol.affine),"para_output_fodf_img.nii" )
with open('fodf_array.npy', 'wb') as f:
    np.save(f, odf)

try:
    shutil.rmtree(path)
except:
    print("Couldn't delete folder")


