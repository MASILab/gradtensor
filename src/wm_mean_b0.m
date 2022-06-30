 r_wm = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/reg_epi/repi_re_fast_wmseg.nii';
rwm_mask = spm_vol(r_wm);
rwm_mask = spm_read_vols(rwm_mask);
rwm_mask(isnan(rwm_mask)) = 0;
wm_mask = logical(rwm_mask);
mean(b0(wm_mask))