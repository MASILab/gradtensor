% simulate the gradient non-linerity in phantom 
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')

%function [simbvec_files,simbval_files,rLimg_file] = sim_gradten_phantom( ...
%	Limg_file, ...
%	refimg_file, ...
%	cor_bval_file, ...
%	cor_bvec_file, ...
%	out_dir ...
%	)

dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.nii';
bvec_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/';
bval_folder = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/';
rLimg_file = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii';
% Load resampled L matrix
VL = spm_vol(rLimg_file);
L = spm_read_vols(VL);
L = reshape(L,[],9);
nv = size(L,1);
vL = zeros(3,3,nv);
vL(1,1,:) = L(:,1);
vL(1,2,:) = L(:,2);
vL(1,3,:) = L(:,3);
vL(2,1,:) = L(:,4);
vL(2,2,:) = L(:,5);
vL(2,3,:) = L(:,6);
vL(3,1,:) = L(:,7);
vL(3,2,:) = L(:,8);
vL(3,3,:) = L(:,9);

% Load data
    %dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
    %b0_vol = dwi_vols(:,:,:,1);
    %dwi_vols = dwi_vols(:,:,:,2:end);

% Get adjbvec and adjbval
    gunzip(fullfile(bvec_folder,'*.gz'))
    bvec_paths_list = dir(fullfile(bvec_folder,'*.nii'));
    bvec_paths = fullfile(bvec_folder, {bvec_paths_list.name});

    bvec_vols = [];
    nb = length(bvec_paths);
    adjbvec = nan(626688,3,nb);
    for i = 1:nb
        load_vol_bvec = spm_vol(bvec_paths{i});
        vol_bvec = spm_read_vols(load_vol_bvec);
        vol_bvec = reshape(vol_bvec(:,:,:,:),[626688,3]);
        adjbvec(:,:,i) = vol_bvec;
        disp(size(vol_bvec));
    end

    gunzip(fullfile(bval_folder,'*.gz'))
    bval_paths_list = dir(fullfile(bval_folder,'*.nii'));
    bval_paths = fullfile(bval_folder, {bval_paths_list.name});
    bval_vols = [];
    adjbval = nan(626688,1,nb);
    for i = 1:nb
        load_vol_bval = spm_vol(bval_paths{i});
        vol_bval = spm_read_vols(load_vol_bval);
        vol_bval = reshape(vol_bval(:,:,:),[626688,1]);
        adjbval(:,i) = vol_bval;
        disp(size(vol_bval));
    end
  
  adjbvec(:,1,:) = -adjbvec(:,1,:);
  bval = nan(1,nb);
  
  ab = adjbvec(1,:,:);
  ab = squeeze(ab);

  %len = sum(ab.^2);
  %len = sum(sqrt(ab));
  %len = len.^2;
  %lenkeeps = len~=0;
  %ab(:,lenkeeps) = ab(:,lenkeeps) .* repmat(len(lenkeeps),3,1);


  len2 = sum(ab.^2);
  %len2 = sum(sqrt(ab));
  bval(1,:) = adjbval(1,:) ./ len2;
  
  bvec = vL(:,:,1).' * ab;

