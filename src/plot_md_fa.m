% plot the uncorrect and correct MD
addpath(genpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home-nfs/masi-shared-home/home/local/VANDERBILT/hansencb/NIFTI');


uncorr_md_posA = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dti_MD.nii.gz');
uncorr_md_posB = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/INPUTS/dti_MD.nii.gz');
uncorr_md_posC = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/INPUTS/dti_MD.nii.gz');


corr_md_posA = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii');
corr_md_posB = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii');
corr_md_posC = load_nii('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/OUTPUTS_future_fieldmap/p_3tb_posC_mask_md.nii');


figure
m = 0 ; 
M = 0.00001 ;
slice = uncorr_md_posA.img(:,47,:);
slice = squeeze(slice);
subplot(1,2,1)
imshow(rot90(slice,1), [m, M], 'Colormap', parula);
title('uncorr md posA');
colorbar;


slice = corr_md_posA.img(:,47,:);
slice = squeeze(slice);
subplot(1,2,2)
imshow(rot90(slice,1), [m, M], 'Colormap', parula);
title('corr md posA');
colorbar;


