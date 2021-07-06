addpath(genpath('/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/matlab/justinlib_v1_7_0'));
addpath('/home/local/VANDERBILT/hansencb/masimatlab/trunk/users/blaberj/dwmri_libraries/dwmri_visualizer_v1_2_0');
addpath('/home/local/VANDERBILT/hansencb/NIFTI');

scan_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/INPUTS';
% scan_dir = '/nfs/masi/hansencb/10_29_2019_human_repositioned/3tb/posB';
mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posC/mask.nii.gz';
dwi_path = fullfile(scan_dir, 'dwmri.nii.gz');
bval_path = fullfile(scan_dir, 'dwmri.bval');
bvec_path = fullfile(scan_dir, 'dwmri.bvec');


dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
mask_vol = logical(nifti_utils.load_untouch_nii4D_vol_scaled(mask_path,'double'));
b0_vol = dwi_vols(:,:,:,1);
dwi_vols = dwi_vols(:,:,:,2:end);

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

est_dwi = zeros(size(dwi_vols));
for i = 1:length(bvals)
    g = bvecs(:,i);
    b = bvals(i);
    
    for x = 1:size(D,1)
        for y = 1:size(D,2)
            for z = 1:size(D,3)
                est_dwi(x,y,z,i) = b0_vol(x,y,z)*exp(-1*b*g'*squeeze(D(x,y,z,:,:))*g);
            end
        end
    end
    
end

diff = abs(est_dwi - dwi_vols);
noise = std(diff,0,4);
snr = nanmean(est_dwi,4) ./ noise;
nanmean(snr(mask_vol))

