addpath(genpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/NIFTI');

none_md_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_MD.nii.gz';
none_fa_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_FA.nii.gz';

est_md_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii';
est_fa_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii';

%mfg_md_path = '/fs4/masi/hansencb/dtiQA_7_test/ICE_20180122_3TA/eddy_drift_sharm_true/501/MD.nii.gz';
%mfg_fa_path = '/fs4/masi/hansencb/dtiQA_7_test/ICE_20180122_3TA/eddy_drift_sharm_true/501/FA.nii.gz';

true_MD = [0.001127; 0.001127; 0.001127; 0.000843; 0.000843; 0.000607; 0.000607; 0.000403; 0.000403; 0.000248; 0.000248; 0.000128; 0.000128];

mask_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii';

mask = nifti_utils.load_untouch_nii_vol(mask_path, 'double');
none_md = nifti_utils.load_untouch_nii_vol(none_md_path, 'double');
none_fa = nifti_utils.load_untouch_nii_vol(none_fa_path, 'double');
est_md = nifti_utils.load_untouch_nii_vol(est_md_path, 'double');
est_fa = nifti_utils.load_untouch_nii_vol(est_fa_path, 'double');

all_vials = mask > 0;
all_vials = logical(all_vials);
none_md(~all_vials) = NaN;
est_md(~all_vials) = NaN;
none_fa(~all_vials) = NaN;
est_fa(~all_vials) = NaN;

for i = 1:13
    vial_mask = logical(mask == i);
    none_md(vial_mask) = none_md(vial_mask) -true_MD(i);
    est_md(vial_mask) = est_md(vial_mask) -true_MD(i);
end

slice = round(size(mask,3)/2);


m = 0;
M = 0.1
figure;
subplot(2,2,1);
imshow(abs(none_fa(:,:,slice)), [m M], 'Colormap', parula);
title('No Correction');
ylabel('FA');
subplot(2,2,2);
imshow(abs(est_fa(:,:,slice)), [m M], 'Colormap', parula);
title('Est. Fields');
%subplot(2,3,3);
%imshow(abs(mfg_fa(:,:,slice)), [m M], 'Colormap', parula);
%title('Mfg. Fields');
colorbar;

m =-0.0001;
m=0;
M = 0.0001
subplot(2,2,3);
imshow(abs(none_md(:,:,slice)), [m M], 'Colormap', parula);
ylabel('MD');
subplot(2,2,4);
imshow(abs(est_md(:,:,slice)), [m M], 'Colormap', parula);
%subplot(2,3,6);
%imshow(abs(mfg_md(:,:,slice)), [m M], 'Colormap', parula);
colorbar;

fa_err_diff_none = abs(none_fa) - abs(est_fa);
%fa_err_diff_mfg = abs(mfg_fa) - abs(est_fa);
md_err_diff_none = abs(none_md) - abs(est_md);
%md_err_diff_mfg = abs(mfg_md) - abs(est_md);


m = -.03;
M = 0.03;
figure;
subplot(1,2,1);
imshow(fa_err_diff_none(:,:,slice), [m M], 'Colormap', parula);
title('Error diff');
ylabel('FA');
%subplot(2,2,2);
%imshow(fa_err_diff_mfg(:,:,slice), [m M], 'Colormap', parula);
%title('Est. Fields');
%title('Mfg. Fields');
colorbar;

m =-0.00005;
M = 0.00005;
subplot(1,2,2);
imshow(md_err_diff_none(:,:,slice), [m M], 'Colormap', parula);
ylabel('MD');
%subplot(2,2,4);
%imshow(md_err_diff_mfg(:,:,slice), [m M], 'Colormap', parula);
colorbar;



