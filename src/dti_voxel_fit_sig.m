% convert bvec and bvec to volumes
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm_vol
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))


% Load data
%dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/rot_45_est_sig.nii';
%dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/posB_3tb_Lest_sig.nii';
dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.nii';
%dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot45_Lest_sig.nii';
dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
b0_vol = dwi_vols(:,:,:,1);
dwi_vols = dwi_vols(:,:,:,2:end);
    
% get the size from dwmri
VD = spm_vol(dwi_path);
D = spm_read_vols(VD);
disp(size(D));
D = reshape(D,[],25);
nv = size(D,1);

% load 
bvec_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bvec';
bval_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bval';
bvec = load(bvec_path);
bvec = bvec();
bval = load(bval_path);
nb = length(bval);
disp(size(bvec))
disp(size(bval))
bval(:,1) = [];
bvec(:,1) = [];

mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii';
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);

% Initialize DT and exitcode volumes
exitcode_vol = zeros(size(b0_vol));
eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
primary_vec_vol = zeros(size(eig_vol)); 
DT_rot_lest_corr_posB = zeros(3,3,size(dwi_vols,1),size(dwi_vols,2),size(dwi_vols,3));
    for i = 1:size(mask_vol,1)
        for j = 1:size(mask_vol,2)
            for k = 1:size(mask_vol,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';
                    %bvecs = squeeze(bvec(i,j,k,:,:))';
                    %bvals = squeeze(bval(i,j,k,:))';

                    % Get linear model
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvec,bval);
                    
                    if sum(isnan(DT_mat(:))) == 0 && sum(isinf(DT_mat(:))) == 0
                        DT_rot_lest_corr_posB(:,:,i,j,k) = DT_mat(:,:);
                    	[v, e] = eig(DT_mat);
                    	e = diag(e);
                    	[max_eig, pos] = max(e);
                    	primary = v(:,pos);
						
                    	eig_vol(i,j,k,:) = e;
						primary_vec_vol(i,j,k,:) = primary;
                        exitcode_vol(i,j,k) = exitcode;
                    end
                end
            end
        end
    end

%save('DT_rot_lest_corr_posB.mat', 'DT_rot_lest_corr_posB');
%out_dir = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/';
%out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_signal_estimate/';
%out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/';
out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
out_name = 'dwmri';
MD = (eig_vol(:,:,:,1) + eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./3;
FA = sqrt(1/2) .* (sqrt( (eig_vol(:,:,:,1) - eig_vol(:,:,:,2)).^2 + (eig_vol(:,:,:,2) - eig_vol(:,:,:,3)).^2  + (eig_vol(:,:,:,3) - eig_vol(:,:,:,1)).^2 ) ./ sqrt(eig_vol(:,:,:,1).^2 + eig_vol(:,:,:,2).^2 + eig_vol(:,:,:,3).^2));
nii = load_untouch_nii(dwi_path);
nii.img = MD;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_md']),nii,'double');
nii.img = FA;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_fa']),nii,'double');
nii.img = primary_vec_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_primary_eigvec']),nii,'double');
   
