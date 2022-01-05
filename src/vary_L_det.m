function [L_det, vL] = vary_L_det()
	% vary_L_det Compute the determinant of LR matrix
	%
	% Outputs
	% L_det		LR determinant
	% vL		LR of the image
% Required packages
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))

% Inputs output dir, mask and resampled LR field
out_dir ='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap';
mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);
rLimg_file = fullfile(out_dir,'L_resamp.nii');

% Reshape the LR field
VL = spm_vol(rLimg_file);
L = spm_read_vols(VL);
size_x = size(L,1);
size_y = size(L,2);
size_z = size(L,3);
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
vL = reshape(vL,[3,3,size_x,size_y,size_z]);
L_det = zeros(size_x,size_y,size_z);

% Along all axis obtain the determinant of LR matrix
for x = 1:size_x
        for y = 1:size_y
            for z = 1:size_z
                %if mask_vol(x,y,z)
                    L_mat = squeeze(vL(:,:,x,y,z));
                    L_det(x,y,z) = det(L_mat(:,:));
                %end
            end
        end
end
