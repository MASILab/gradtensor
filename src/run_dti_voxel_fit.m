% run dti_voxel_fit on the dwmri and corrected bvec and bvals 

% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')

%gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii.gz')

dwi_path = '/home/local/VANDERBILT/kanakap/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.nii';
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/corrected_bval/';
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
%mask_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwi_mask.nii';
bvec_file = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bvec';
bval_file = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bval';
mask_path = '';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_real_masiver/full/';
out_name = 'real_corr_bx';

%dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
dti_voxel_fit_masiver(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name,bval_file,bvec_file);
