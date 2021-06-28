import numpy as np
import nibabel as nib
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.reconst.dti import (TensorModel, color_fa, fractional_anisotropy,
                              geodesic_anisotropy, mean_diffusivity,
                              axial_diffusivity, norm,
                              radial_diffusivity, lower_triangular)

in_dwi = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.nii.gz'
img = nib.load(in_dwi)
data = img.get_fdata(dtype=np.float32)
affine = img.affine

in_bval = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.bval'
in_bvec = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.bvec'

outfile = '/home/local/VANDERBILT/kanakap/INPUTS/recon_norm_dwmri.nii.gz'
diff_file = '/home/local/VANDERBILT/kanakap/INPUTS/recon_diff_norm_dwmri.nii.gz'
p_diff_file = '/home/local/VANDERBILT/kanakap/INPUTS/recon_pdiff_norm_dwmri.nii.gz'

bvals, bvecs = read_bvals_bvecs(in_bval, in_bvec)
gtab = gradient_table(bvals, bvecs, b0_threshold=bvals.min())

#in_mask = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii.gz'
#mask = get_data_as_mask(nib.load(in_mask), dtype=bool)

def _get_min_nonzero_signal(data):
    return np.min(data[data > 0])

S0 = np.mean(data[..., gtab.b0s_mask], axis=-1)
B0 = data[:,:,:,0]
diffus_sig = np.zeros(data.shape, dtype=np.float32)
data_diff = np.zeros(data.shape, dtype=np.float32)
data_diff_percent = np.zeros(data.shape, dtype=np.float32)

norm_data = np.zeros(data.shape, dtype=np.float32)

tenmodel = TensorModel(gtab, fit_method='LS',
                               min_signal=_get_min_nonzero_signal(data))

def noramlize_diff():
    for i in range(data.shape[3]):
        norm_data[:,:,:,i] = data[:,:,:,i] / B0[:,:,:]
    return norm_data


for i in range(data.shape[0]):
   # if mask is not None:
    #    tenfit2 = tenmodel.fit(data[i, :, :, :], mask[i, :, :])
    #else:
    tenfit2 = tenmodel.fit(data[i, :, :, :])
    print(i)
        #data_diff[i, :, :, :] = np.abs(tenfit2.predict(gtab, S0[i, :, :]).astype(np.float32) - data[i, :, :])
    diffus_sig[i, :, :, :] = np.abs(tenfit2.predict(gtab, S0[i, :, :]).astype(np.float32))
    DS = diffus_sig[..., ~gtab.b0s_mask]
    DS_img = nib.Nifti1Image(DS.astype(np.float32), affine)
    nib.save(DS_img, outfile)

    norm_data = noramlize_diff()
    data_diff[i, :, :, :] = np.abs(tenfit2.predict(gtab, S0[i, :, :]).astype(np.float32) - norm_data[i, :, :])
    R = data_diff[..., ~gtab.b0s_mask]
    R_img = nib.Nifti1Image(R.astype(np.float32), affine)
    nib.save(R_img, diff_file)
    
    data_diff_percent[i, :, :, :] = 100*(np.abs(tenfit2.predict(gtab, S0[i, :, :]).astype(np.float32) - norm_data[i, :, :]))/norm_data[i, :, :]
    PR = data_diff_percent[..., ~gtab.b0s_mask]
    PR_img = nib.Nifti1Image(PR.astype(np.float32), affine)
    nib.save(PR_img, p_diff_file)
