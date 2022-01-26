
Test_dwi(:,:,:,i) = abs(est_dwi(:,:,:,i) +  real_noise(i) + img_noise(i));

for i = 5:24
    Test_dwi(:,:,:,i) = abs(est_dwi(:,:,:,i) +  real_noise(i) + img_noise(i));
end

mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii' ;
 mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);


est = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_realdata100/est_fa.nii';
Lest = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_realdata100/Lest_fa.nii';
NLest = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_realdata100/NLest_fa.nii';
true = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii';
est_v = nifti_utils.load_untouch_nii4D_vol_scaled(est,'double');
Lest_v = nifti_utils.load_untouch_nii4D_vol_scaled(Lest,'double');
NLest_v = nifti_utils.load_untouch_nii4D_vol_scaled(NLest,'double');
true = nifti_utils.load_untouch_nii4D_vol_scaled(true,'double');

 diff = est_v - true;
Ldiff = Lest_v - true;
nanmean(nanmean(nanmean(diff(mask_vol))));
nanmean(nanmean(nanmean(Ldiff(mask_vol))));
 NLdiff = NLest_v - true;
rNLdiff = reshape(NLdiff, [], 1);
%rNLest_v = reshape(NLest_v, [], 1);
nanmean(nanmean(nanmean(NLest_v))) / std(rNLdiff, 'omitnan');


%% final testing
NLest = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_realdata/NLest_fa.nii';
true = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii';
NLest_v = nifti_utils.load_untouch_nii4D_vol_scaled(NLest,'double');
true = nifti_utils.load_untouch_nii4D_vol_scaled(true,'double');
 NLdiff = NLest_v - true;
rNLdiff = reshape(NLdiff, [], 1);
%rNLest_v = reshape(NLest_v, [], 1);
nanmean(nanmean(nanmean(NLest_v))) / std(rNLdiff, 'omitnan');
