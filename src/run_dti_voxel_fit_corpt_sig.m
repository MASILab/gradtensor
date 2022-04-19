% run dti_voxel_fit_sig on the corpt dwmri and original bvec bvals 

% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm_vol
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))

%gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii.gz')

dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/BLSA_Lest_sig.nii';
%dwi_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwisubj2scanner.nii'
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
mask_path = '';
bvec_path = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bvec';
bval_path = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bval'
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver/';
%out_dir = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/'
%out_name = 'Lest';
out_name = 'masiver'

%dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
dti_voxel_fit_sig(dwi_path,bvec_path,bval_path,mask_path, out_dir, out_name);
