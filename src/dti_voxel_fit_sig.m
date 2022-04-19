function dti_voxel_fit_sig(dwi_path,bvec_path,bval_path,mask_path, out_dir, out_name)
    % dti_voxel_fit_sig Compute FA, MD, PEV for approximated corrected signal
    %
    % Inputs
    %   dwi_path        Path to the diffusion weighted nifti. Assumes b0 at first volume
    %   bvec_path      Cell array of paths to bvec niftis. Assumes these are with respect to voxels.
    %   bval_path      Ccell array of paths to bval niftis.
    %   mask_path       Optional path to mask nifti.
    %   out_dir         Path to a directory in which to save metrics
    %   out_name        Prefix to the generated metric nifti filenames

% Load the required packages
% for nifti_utils
%addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
%addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm_vol
%addpath(genpath('../external/spm_read_nii'))
%addpath(genpath('../external/spm_reslice'))

% Load dwmri data
dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
bvec = load(bvec_path);
bval = load(bval_path);

% dwi and non-dwi
ind_b0 = find(~bval);
ind_non_b0 = find(bval);

all_b0_vol = dwi_vols(:,:,:,ind_b0);
b0_vol = mean(all_b0_vol,4);
dwi_vols = dwi_vols(:,:,:,ind_non_b0);
    
% Load bvec and bval 
bvec = load(bvec_path);
bvec = bvec();
bval = load(bval_path);
nb = length(bval);
disp(size(bvec))
disp(size(bval))
bval = bval(ind_non_b0);
bvec = bvec(:,ind_non_b0);

% Load mask 
%mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
if ~exist('mask_path','var') || isempty(mask_path)
        disp('making mask')
        mask_vol = true(size(b0_vol));
    else
        %gunzip(mask_path)
        mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
        mask_vol = logical(mask_vol);
    end
%mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
%mask_vol = logical(mask_vol);

% Initialize DT and exitcode volumes
exitcode_vol = zeros(size(b0_vol));
eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
primary_vec_vol = zeros(size(eig_vol)); 
AD = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3));
RD = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3));

% Cycle over and compute DT voxel-wise
    for i = 1:size(dwi_vols,1)
        for j = 1:size(dwi_vols,2)
            for k = 1:size(dwi_vols,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';

                    % Get linear model
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvec,bval);
                    
                    if sum(isnan(DT_mat(:))) == 0 && sum(isinf(DT_mat(:))) == 0
                        %DT_rot_lest_corr_posB(:,:,i,j,k) = DT_mat(:,:);
                    	[v, e] = eig(DT_mat);
                    	e = diag(e);
                    	[max_eig, pos] = max(e);
                    	primary = v(:,pos);

                        [min_eig, ter_pos] = min(e);
                        possible_positions = [1 2 3];
                        not_max = possible_positions(possible_positions~=pos);
                        sec_pos = not_max(not_max~=ter_pos);
                        sec_eig = e(sec_pos);
						
                    	eig_vol(i,j,k,:) = e;
                        primary_vec_vol(i,j,k,:) = primary;

                        AD(i,j,k) = max_eig;
                        RD(i,j,k) = (min_eig + sec_eig) / 2;

                        exitcode_vol(i,j,k) = exitcode;
                    end
                end
            end
        end
    end

% Compute FA, MD, PEV and save it
%out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_noise/';
%out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_noise100/';
%out_name = 'Lest_corr_sm';	
%out_name = 'Lest'
MD = (eig_vol(:,:,:,1) + eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./3;
FA = sqrt(1/2) .* (sqrt( (eig_vol(:,:,:,1) - eig_vol(:,:,:,2)).^2 + (eig_vol(:,:,:,2) - eig_vol(:,:,:,3)).^2  + (eig_vol(:,:,:,3) - eig_vol(:,:,:,1)).^2 ) ./ sqrt(eig_vol(:,:,:,1).^2 + eig_vol(:,:,:,2).^2 + eig_vol(:,:,:,3).^2));
%AD = eig_vol(:,:,:,1);
%RD = (eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./2;
nii = load_untouch_nii(dwi_path);
nii.img = MD;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_md']),nii,'double');
nii.img = FA;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_fa']),nii,'double');
nii.img = primary_vec_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_primary_eigvec']),nii,'double');
nii.img = AD;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_ad']),nii,'double');
nii.img = RD;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_rd']),nii,'double');   
end
