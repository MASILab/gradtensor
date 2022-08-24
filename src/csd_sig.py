import os
import nibabel as nib
import numpy as np
import subprocess

from dipy.reconst.csdeconv import ConstrainedSphericalDeconvModel
from dipy.data import get_sphere
from dipy.core.gradients import gradient_table
from scilpy.reconst.multi_processes import fit_from_model, convert_sh_basis

n = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/sub-cIVs001_ses-s1Bx2_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.nii').get_fdata()
vol = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/sub-cIVs001_ses-s1Bx2_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.nii')
mask = nib.load('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/mask.nii.gz').get_fdata()
vec_folder = '/home/local/VANDERBILT/kanakap/try_fx_emp/bvec'
val_folder = '/home/local/VANDERBILT/kanakap/try_fx_emp/bval'

bvec_vols = []
for i in sorted(os.listdir(vec_folder)):
    if i.endswith('.nii'):
        bvec_vol = nib.load(vec_folder + '/' + i).get_fdata()
        bvec_vol = np.expand_dims(bvec_vol,4)
        bvec_vol = np.transpose(bvec_vol,(0,1,2,4,3))
        bvec_vols.append(bvec_vol)
bvec_stack = np.stack(bvec_vols,3)
bvec_stack = bvec_stack.squeeze()

bval_vols = []
for i in sorted(os.listdir(val_folder)):
    if i.endswith('.nii'):
        bval_vol = nib.load(val_folder + '/' + i).get_fdata()
        bval_vols.append(bval_vol)
bval_stack = np.stack(bval_vols,3)

reg_sphere = get_sphere('symmetric362')
odf = np.zeros((n.shape[1],n.shape[2],n.shape[3],45))
for i in range(112):
    for j in range(112):
        for k in range(54):
            if mask[i,j,k] == 1:
                vec = bvec_stack[i,j,k,:,:]
                val = bval_stack[i,j,k,:]
                dwi = n[i][j][k]
                dwi = np.expand_dims(dwi,(1,2,3))
                dwi = np.transpose(dwi,(1,2,3,0))
                gtab = gradient_table(val, vec)
                csd_model = ConstrainedSphericalDeconvModel(gtab, ((1.509867245705930885e-03,3.422135009646922127e-04,3.422135009646922127e-04), 1.058279199218750000e+04),reg_sphere=reg_sphere,sh_order=8)
                csd_fit = fit_from_model(csd_model, dwi)
                shm_coeff = csd_fit.shm_coeff
                shm_coeff = convert_sh_basis(shm_coeff, reg_sphere)
                odf[i,j,k,:] = shm_coeff
nib.save(nib.Nifti1Image(odf.astype(np.float32),vol.affine),"output_fodf_img.nii" )

            vec = np.array(vec)
            tvec = np.transpose(vec)
            val = np.array(val)
            tval = np.reshape(val,[1,112])
            dwi.append(n[i][j][k])
            
            bashCommand = "scil_compute_ssst_frf.py dwi tval tvec"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

#csd_model = ConstrainedSphericalDeconvModel(gtab, ((1.509867245705930885e-03,3.422135009646922127e-04,3.422135009646922127e-04), 1.058279199218750000e+04),reg_sphere=reg_sphere,sh_order=8)           
# tired
# from scilpy.reconst.frf import compute_ssst_frf
# resp = compute_ssst_frf(dwi,val,vec)
