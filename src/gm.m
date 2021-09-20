% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm_vol
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))


% Load data
gunzip( '/nfs/masi/hansencb/10_29_2019_human_repositioned/3ta/posB/slant/OUTPUTS/FinalResult/T1_seg.nii.gz');
seg_path ='/nfs/masi/hansencb/10_29_2019_human_repositioned/3ta/posB/slant/OUTPUTS/FinalResult/T1_seg.nii';
%dwi_vols = nifti_utils.load_untouch_nii(seg_path,'double');
dwi_vols = load_untouch_nii(seg_path);
dwi_vols = dwi_vols.img;
%dwi_vols = niftiread(seg_path)
dwi_vols(dwi_vols==40)  = 0;
dwi_vols(dwi_vols==41)  = 0;
dwi_vols(dwi_vols==44) = 0;
dwi_vols(dwi_vols==45) = 0;
dwi_vols(dwi_vols==4) = 0;
dwi_vols(dwi_vols==11) = 0;
dwi_vols(dwi_vols==49) = 0;
dwi_vols(dwi_vols==50) = 0;
dwi_vols(dwi_vols==51) = 0;
dwi_vols(dwi_vols==52) = 0;
%dwi_vols(dwi_vols==40 & dwi_vols==41 & dwi_vols==44 & dwi_vols==45 & dwi_vols==4 & dwi_vols==11 & dwi_vols==49 & dwi_vols==50  & dwi_vols==51 & dwi_vols==52) = 0;
dwi_vols(dwi_vols>=1) = 1;

out_name = 'slant';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB';
%nii = niftiwrite(dwi_vols, '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/slant_gm_seg1.nii');
nii = load_untouch_nii(seg_path);
nii.img = dwi_vols;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_gm_seg2']),nii,'double');
