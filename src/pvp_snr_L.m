addpath(genpath('/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home/local/VANDERBILT/hansencb/NIFTI');

scan_dir = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/';
% scan_dir = '/nfs/masi/hansencb/10_29_2019_human_repositioned/3tb/posB';
mask_path = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii.gz';
out_dir  = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/simulate_brain/';
dwi_path = fullfile(scan_dir, 'dwmri.nii.gz');
bval_path = fullfile(scan_dir, 'dwmri.bval');
bvec_path = fullfile(scan_dir, 'dwmri.bvec');

dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
mask_vol = logical(nifti_utils.load_untouch_nii4D_vol_scaled(mask_path,'double'));
b0_vol = dwmri_vols(:,:,:,1);
dwi_vols = dwmri_vols(:,:,:,2:end);

bvecs = importdata(bvec_path);
bvals = importdata(bval_path);
bvecs = bvecs(:,2:end);
bvals = bvals(2:end);

[DT_vol, tmp] = linear_vol_fit(b0_vol, dwi_vols, bvecs, bvals, mask_vol);

D = zeros(size(DT_vol,1),size(DT_vol,2),size(DT_vol,3),3,3);
D(:,:,:,1,1) = DT_vol(:,:,:,1);
D(:,:,:,1,2) = DT_vol(:,:,:,2);
D(:,:,:,1,3) = DT_vol(:,:,:,3);
D(:,:,:,2,1) = DT_vol(:,:,:,2);
D(:,:,:,2,2) = DT_vol(:,:,:,4);
D(:,:,:,2,3) = DT_vol(:,:,:,5);
D(:,:,:,3,1) = DT_vol(:,:,:,3);
D(:,:,:,3,2) = DT_vol(:,:,:,5);
D(:,:,:,3,3) = DT_vol(:,:,:,6);

gunzip('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii.gz')
rL_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii';
% Load resampled L
VL = spm_vol(rL_path);
L = spm_read_vols(VL);
%L = reshape(L,[],9);
nv = size(L,1);
vL = zeros(3,3,size(L,1),size(L,2),size(L,3));
vL(1,1,:,:,:) = L(:,:,:,1);
vL(1,2,:,:,:) = L(:,:,:,2);
vL(1,3,:,:,:) = L(:,:,:,3);
vL(2,1,:,:,:) = L(:,:,:,4);
vL(2,2,:,:,:) = L(:,:,:,5);
vL(2,3,:,:,:) = L(:,:,:,6);
vL(3,1,:,:,:) = L(:,:,:,7);
vL(3,2,:,:,:) = L(:,:,:,8);
vL(3,3,:,:,:) = L(:,:,:,9);

est_dwi = zeros(size(dwi_vols));
for i = 1:length(bvals)
    g = bvecs(:,i);
    b = bvals(i);
    
    for x = 1:size(D,1)
        for y = 1:size(D,2)
            for z = 1:size(D,3)
                est_dwi(x,y,z,i) = b0_vol(x,y,z)*exp(-1*b*g'*squeeze(vL(:,:,x,y,z))'*squeeze(D(x,y,z,:,:))*squeeze(vL(:,:,x,y,z))*g);
            end
        end
    end    
end

out_name = 'org';
dwmri_est = zeros(size(dwmri_vols));
dwmri_est(:,:,:,1) = b0_vol ;
dwmri_est(:,:,:,2:end) = est_dwi ; 
nii = load_untouch_nii(dwi_path);
nii.img = dwmri_est;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_sig']),nii,'double');

diff = abs(est_dwi - dwi_vols);
noise = std(diff,0,4);
snr = nanmean(est_dwi,4) ./ noise;
nanmean(snr(mask_vol))
