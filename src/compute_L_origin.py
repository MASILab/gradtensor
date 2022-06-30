import nibabel as nib
import numpy as np
import math
img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS_zero_grad/L.nii.gz')
L = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS_zero_grad/L.nii.gz').get_fdata()
OL = np.where(L==math.nan,L,1)
OL_img = nib.Nifti1Image(OL, img.affine, img.header)
nib.save(OL_img, '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS_zero_grad/L_origin.nii')



img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/B0_X_Hz.nii')
 ones = np.ones([96,96,96])
x = nib.Nifti1Image(ones, img.affine, img.header)
nib.save(x, '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/X_linear.nii')

img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/B0_Y_Hz.nii')
y = nib.Nifti1Image(ones, img.affine, img.header)
nib.save(y, '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/Y_linear.nii')

img = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/B0_Z_Hz.nii')
z = nib.Nifti1Image(ones, img.affine, img.header)
nib.save(z, '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_true_fieldmap/INPUTS/Z_linear.nii')


