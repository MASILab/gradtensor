addpath(genpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/NIFTI');


est_md_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii';
mask_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii';

mask = nifti_utils.load_untouch_nii_vol(mask_path, 'double');
est_md = nifti_utils.load_untouch_nii_vol(est_md_path, 'double');


slice = round(size(mask,3)/2);


%m = -0.0001;
m = 0;
M = 0.0001;

figure;
%subplot(2,3,1);

imshow(abs(est_md(:,:,slice)), [m M], 'Colormap', parula);
title('Est Correction');
ylabel('MD');
colorbar;