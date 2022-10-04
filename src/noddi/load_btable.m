
% addpath '/Users/kanakap/VUFall/DWI Project/code/HowToRenderVideo-master/matlab'
addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');
bvec_folder = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_lr_corr_1/Lemp/emp_Lcorrected_bvec';
bval_folder = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_lr_corr_1/Lemp/emp_Lcorrected_bval';
bval_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/sub-cIVs001_ses-s1Bx2_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.bval';

% Remove the b0 vols
bval = load(bval_file);
non_b0_idx = find(bval ~= 0);

% Exctract nii and create cell array of bvec nifti paths
gunzip(fullfile(bvec_folder,'*.gz'))
bvec_paths_list = dir(fullfile(bvec_folder,'*.nii'));
bvec_paths = fullfile(bvec_folder, {bvec_paths_list.name});    

bvec_vols = [];
for i = 1:length(bvec_paths)
    bvec_vols = cat(4,bvec_vols,permute(nifti_utils.load_untouch_nii_vol_scaled(bvec_paths{i},'double'),[1 2 3 5 4]));
end
% bvec_vols = bvec_vols(:,:,:,non_b0_idx,:);

% Exctract nii and create cell array of bval nifti paths
gunzip(fullfile(bval_folder,'*.gz'))
bval_paths_list = dir(fullfile(bval_folder,'*.nii'));
bval_paths = fullfile(bval_folder, {bval_paths_list.name});

bval_vols = [];
for i = 1:length(bval_paths)
    bval_vols = cat(4,bval_vols,nifti_utils.load_untouch_nii_vol_scaled(bval_paths{i},'double'));
end
% bval_vols = bval_vols(:,:,:,non_b0_idx);

