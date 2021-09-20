% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))
flags = struct( ...
        'mask',true, ...
        'mean',false, ...
        'interp',1, ...
        'which',1, ...
        'wrap',[0 0 0], ...
        'prefix','r' ...
        );
refimg_file = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.nii'
gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/slant_wm_seg.nii.gz')
Limg_file = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/slant_wm_seg.nii'
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB'
spm_reslice({refimg_file; Limg_file},flags);
[p,n,e] = fileparts(Limg_file);
rLimg_file = fullfile(p,['r' n e]);
movefile(rLimg_file,fullfile(out_dir,'wm_seg_resamp.nii'));
rLimg_file = fullfile(out_dir,'wm_seg_resamp.nii');
