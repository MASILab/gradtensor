% run compute_corrput_signal - to compute the dwi signal with and without corrupt

% required packages
% for nifti_utils
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))


our_dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.nii';
gunzip('/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwisubj2scanner.nii.gz')
%dwi_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwisubj2scanner.nii';
dwi_path = '/home/local/VANDERBILT/kanakap/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.nii';
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/';
%org_bvec_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bvec';
%org_bval_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bval';
org_bvec_path = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bvec';
org_bval_path = '/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIIIs01/ses-s1Bx1/dwi/sub-cIIIs01_ses-s1Bx1_acq-b1000n40r21x21x22peAPP_run-105_dwi.bval';
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
%gunzip('/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwi_mask.nii.gz')
mask_path = '';
%mask_path = '/nfs/masi/kanakap/projects/LR/population_basis_study/data/reg/dwi_mask.nii';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_masiver';
out_name = 'BLSA';
%rL_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii';
L_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz'
compute_corrput_signal(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name, L_path,org_bvec_path,org_bval_path)

