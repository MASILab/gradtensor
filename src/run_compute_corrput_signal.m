% run compute_corrput_signal - to compute the dwi signal with and without corrupt

% required packages
% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))


%dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.nii';
dwi_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/Diffusion.nii'
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/';
%org_bvec_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bvec';
%org_bval_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bval';
org_bvec_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/Diffusion.bvecs'
org_bval_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/Diffusion.bvals'
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
gunzip('/nfs/masi/kanakap/projects/LR/population_basis_study/data/mask.nii.gz')
mask_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/mask.nii'
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/BLSA_estimates';
out_name = 'BLSA';
rL_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii';

compute_corrput_signal(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name, rL_path,org_bvec_path,org_bval_path)
