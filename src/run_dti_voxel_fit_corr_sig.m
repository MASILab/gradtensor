% run dti_voxel_fit on the approx corrected dwmri and orginal bvec and bvals 

% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm_vol
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))

%gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii.gz')

%dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/approx/real_sig.nii';
dwi_path = '/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.nii';
%bvec_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bvec';
%bval_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bval';
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
bvec_path = '/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bvec';
bval_path = '/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bval';
mask_path = '';
%out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver';
out_dir = '/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_'
%out_name = 'real_corr_sm';
out_name = 'true_'

%dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
dti_voxel_fit_sig(dwi_path,bvec_path,bval_path,mask_path, out_dir, out_name);
