% run dti_voxel_fit on the dwmri and corrected bvec and bvals 

% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')

gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii.gz')

dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/BLSA_estimates/BLSA_Lest_sig.nii';
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/BLSA_bx/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/BLSA_bx/corrected_bval/';
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
mask_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwi_mask.nii';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/BLSA_bx/';
out_name = 'Lest_corr_bx';

%dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
