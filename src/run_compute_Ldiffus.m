% run compute_L_diffus_signal - to compute the dwi signal with and without corrupt

% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))

%gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii.gz')

dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.nii';
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/corrected_bval/';
org_bvec_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bvec';
org_bval_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bval';
mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_estimates_noflip';
out_name = 'posB';
rL_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/L_resamp.nii';

%gunzip(mask_path)
%gunzip(rL_path)
%dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
compute_L_diffus_sig(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name, rL_path,org_bvec_path,org_bval_path)
