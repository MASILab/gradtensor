import nibabel as nib
import numpy as np
img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/slant_T1_seg.nii.gz')
slant_np = img.get_fdata()
slant_np = np.where(slant_np ==40.0, 0.0,slant_np)
slant_np = np.where(slant_np ==41.0, 0.0,slant_np)
slant_np = np.where(slant_np ==44.0, 0.0,slant_np)
slant_np = np.where(slant_np ==45.0, 0.0,slant_np)
slant_np = np.where(slant_np ==4.0, 0.0,slant_np)
slant_np = np.where(slant_np ==11.0, 0.0,slant_np)
slant_np = np.where(slant_np ==49.0, 0.0,slant_np)
slant_np = np.where(slant_np ==50.0, 0.0,slant_np)
slant_np = np.where(slant_np ==51.0, 0.0,slant_np)
slant_np = np.where(slant_np ==52.0, 0.0,slant_np)
slant_np = np.where(slant_np > 0.0 , 1.0,slant_np)
nib.save(nib.Nifti1Image(slant_np,img.affine),'/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/slant_gm_seg.nii.gz')

