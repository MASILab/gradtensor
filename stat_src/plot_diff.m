% plot only MD difference
addpath(genpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/NIFTI');

img_path1 = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_MD.nii.gz';
img_path2 = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii';
mask_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii';
save_path = '/home/local/VANDERBILT/kanakap/INPUTS/diff_posA_MD.nii.gz';

nii1 = load_nii(img_path1);
nii2 = load_nii(img_path2);
mask= load_nii(mask_path);
  
diff = nii1.img - nii2.img;
diff = abs(diff);
nii1.img = diff .* logical(mask.img);
nii1.img(~logical(mask.img)) = 0;
save_nii(nii1, save_path);

%m = -0.0001;
m = 0;
M = 0.0001;
nii3 = load_nii(save_path);
mask = logical(mask.img);
nii3.img(~mask) = m;
figure;
slice = nii3.img(:,round(size(nii3.img,3)/2),:);
slice = squeeze(slice);
imshow(rot90(slice,1), [m, M], 'Colormap', parula);
title('MD Diff');
colorbar;

