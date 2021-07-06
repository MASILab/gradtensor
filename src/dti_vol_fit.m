%% Set environment
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')


dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
bvec_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
bval_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/';

%% Fit diffusion tensor
dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(dwi_path,'dwmri.nii'),'double');
bvecs = dlmread(fullfile(bvec_path,'dwmri.bvec'));
bvals = dlmread(fullfile(bval_path,'dwmri.bval'));
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(mask_path,'mask.nii'),'logical');
 
% Split into B0 and DWI
b0_vol = dwmri_vol(:,:,:,bvals == 0);
dwi_vol = dwmri_vol(:,:,:,bvals ~= 0);
bvecs_dwi = bvecs(:,bvals ~= 0);
bvals_dwi = bvals(bvals ~= 0);
 
% Make output directory
output_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
 
% Loop over and compute diffusion tensor, V1, and FA
V1_vol = zeros(size(dwmri_vol,1),size(dwmri_vol,2),size(dwmri_vol,3),3);
%FA_vol = zeros(size(dwmri_vol,1),size(dwmri_vol,2),size(dwmri_vol,3));
FA_vol = zeros(dwi_vol);
MD_vol = zeros(size(dwmri_vol,1),size(dwmri_vol,2),size(dwmri_vol,3));
for i = 1:size(dwmri_vol,1)
    for j = 1:size(dwmri_vol,2)
        for k = 1:size(dwmri_vol,3)
            if mask_vol(i,j,k)
                % Compute A: -b*g*g'
                gx = bvecs_dwi(1,:)';
                gy = bvecs_dwi(2,:)';
                gz = bvecs_dwi(3,:)';
                A = -bsxfun(@times,bvals_dwi',[gx.^2 2*gx.*gy 2*gx.*gz gy.^2 2*gy.*gz gz.^2]);
 
                % Fit using linear least squares
                DT = mldivide(A,squeeze(log(abs(dwi_vol(i,j,k,:)./b0_vol(i,j,k)))));
                DT_mat = zeros(3,3);
 
                % Reshape to matrix
                DT_mat(1,1:3) = DT(1:3);
                DT_mat(2,2:3) = DT(4:5);
                DT_mat(3,3) = DT(6);
                DT_mat(2,1) = DT_mat(1,2);
                DT_mat(3,1:2) = DT_mat(1:2,3)';
 
                if all(isfinite(DT_mat(:)))
                    % Get eigenvectors and eigenvalues
                    [V,D] = eig(DT_mat);

		    MA_vol(i,j,k) = (D(1) + D(2) + D(3))/3;
                    % Sort by size
                    [D,I] = sort(diag(D),'descend');
                    V = V(:, I);
 
                    % Store V1
                    V1_vol(i,j,k,:) = V(:,1);
 
                    % Compute FA
                    FA_vol(i,j,k,:) = sqrt(1/2)*sqrt((D(1)-D(2))^2 + (D(2)-D(3))^2 + (D(3)-D(1))^2)/sqrt(D(1)^2+D(2)^2+D(3)^2);
		    %MA_vol(i,j,k) = (D(1) + D(2) + D(3))/3;
                end
            end
        end
    end
end 
 
%% Save output
template_nii = load_untouch_nii(fullfile(dwi_path,'dwmri.nii'));
 
% V1
template_nii.img = V1_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(output_dir,'dti_V1'),template_nii,'double');
 
% FA
template_nii.img = FA_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(output_dir,'dti_FA'),template_nii,'double');

%MD
template_nii.img = MD_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(output_dir,'dti_MD'),template_nii,'double');

