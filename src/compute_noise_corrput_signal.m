function compute_noise_corrput_signal(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name, rL_path, org_bvec_path, org_bval_path, initial_SNR)
    % compute_noise_corrput_signal Computes signal with LR and noise induced in it 
    %
    % Inputs
    %   dwi_path 	Path to the diffusion weighted nifti. Assumes b0 at first volume
    %   bvec_folder 	Folder with the corrected voxel wise bvec files 
    %   bval_folder 	Folder with the corrected voxel wise bval files
    %   mask_path 	Path to mask nifti (optional)
    %   out_dir 	Path to a directory in which to save metrics
    %   out_name 	Prefix to the generated metric nifti filenames
    %   rL_path 	Path to resampled L
    %   org_bvec_path 	Orignial dwmri bvec
    % 	org_bvec_path 	Original dwmri bval
    
    % Load data
    dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
    b0_vol = dwmri_vols(:,:,:,1);
    dwi_vols = dwmri_vols(:,:,:,2:end);
    % Load resampled L
    VL = spm_vol(rL_path);
    L = spm_read_vols(VL);
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
    %save('L_mat','vL');

    % Exctract nii and create cell array of bvec nifti paths
    gunzip(fullfile(bvec_folder,'*.gz'))
    bvec_paths_list = dir(fullfile(bvec_folder,'*.nii'));
    bvec_paths = fullfile(bvec_folder, {bvec_paths_list.name});
    
    bvec_vols = [];
    for i = 1:length(bvec_paths)
        bvec_vols = cat(4,bvec_vols,permute(nifti_utils.load_untouch_nii_vol_scaled(bvec_paths{i},'double'),[1 2 3 5 4]));
    end
    
    % Exctract nii and create cell array of bval nifti paths
    gunzip(fullfile(bval_folder,'*.gz'))
    bval_paths_list = dir(fullfile(bval_folder,'*.nii'));
    bval_paths = fullfile(bval_folder, {bval_paths_list.name});
    
    bval_vols = [];
    for i = 1:length(bval_paths)
        bval_vols = cat(4,bval_vols,nifti_utils.load_untouch_nii_vol_scaled(bval_paths{i},'double'));
    end

    % Load the mask nii file
    if ~exist('mask_path','var') || isempty(mask_path)
        mask_vol = true(size(b0_vol));
    else
        %gunzip(mask_path)
        mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
        mask_vol = logical(mask_vol);
    end
    
    % Volume 2 to end for bvec and bval
    bvec_vols = bvec_vols(:,:,:,2:end,:);
    bval_vols = bval_vols(:,:,:,2:end);
    
    % Load the original bvec and bval 
    org_bvecs = importdata(org_bvec_path);
    org_bvals = importdata(org_bval_path);
    org_bvecs = org_bvecs(:,2:end);
    org_bvals = org_bvals(2:end);

    % Initialize DT and exitcode volumes
    exitcode_vol = zeros(size(b0_vol));
    eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
    primary_vec_vol = zeros(size(eig_vol)); 
    est_dwi = zeros(size(dwi_vols));
    Lest_dwi = zeros(size(dwi_vols));
    Nest_dwi = zeros(size(dwi_vols));
    NLest_dwi = zeros(size(dwi_vols));
    
    % WM in b0
    % wm = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/reg_epi/epi_re_fast_wmseg.nii
    %spm_reslice({dwi_path;wm},flags)
    r_wm = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/reg_epi/repi_re_fast_wmseg.nii';
    rwm_mask = spm_vol(r_wm);
    rwm_mask = spm_read_vols(rwm_mask);
    rwm_mask(isnan(rwm_mask)) = 0;
    wm_mask = logical(rwm_mask);

    % Add complex guassian noise
    SNR = initial_SNR;
    mean_signal = mean(b0_vol(wm_mask))
    noise_std = mean_signal / initial_SNR;
    randn('state',sum(100*clock));
    real_noise = noise_std * randn(size(dwi_vols));
    img_noise = 1i * noise_std * randn(size(dwi_vols));
    
    % Create complex guassian noise
    for i = 1:size(mask_vol,1)
        for j = 1:size(mask_vol,2)
            for k = 1:size(mask_vol,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals   
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';
                    bvecs = squeeze(bvec_vols(i,j,k,:,:))';
                    bvals = squeeze(bval_vols(i,j,k,:))';
                    L_mat = squeeze(vL(:,:,i,j,k));

                    % Get linear model - use the corrected bvec and bval - considered as ground 
                    % truth tensor
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvecs,bvals);
                   
                    % Induce corrpution (script adpated from compute_b_images.m) and signal equation
                    for v = 1:length(bvals)
                        g = bvecs(:,v);
                        b = bvals(v);
                        
                        og = org_bvecs(:,v);
                        ob = org_bvals(v);
                        
                       
                        % DWI signal with no corput
                        est_dwi(i,j,k,v) =  b0_vol(i,j,k)*exp(-1*ob*og'*DT_mat(:,:)*og);
           		Nest_dwi(i,j,k,v) = abs(est_dwi(i,j,k,v) + real_noise(i,j,k,v) + img_noise(i,j,k,v)) ;
                 
                        % To compute ADC
                        %bv_b0 = 0;
                        %ADC = (log(est_dwi(i,j,k,v) / (b0_vol(i,j,k))) * (1 / (bv_b0 - ob)));
                        %fprintf('ADC no corpt %f for volume %i\n',[ADC, v]);
                        
                        % Here b and g need to be original uncorrected bvals, and bvecs
                        % Most simply, the adjusted bvec is simply L * bvec. Here we are
                        % operating in the image space.
                        og(1) = -og(1);
                        ab = L_mat * og;

                        % The bvecs were length 1 before adjustment, so now compute the length
                        % change and adjust bvals accordingly. NOTE: adjust bval by the square
                        % of the length, because the b value has a G^2 term but the vector
                        % length is for G.
                        len2 = sum(ab.^2);
                        adjbval = ob .* len2;

                        % Re-normalize bvecs to length 1 to compensate for the b value
                        % adjustment we just made. Skip cases where b=0.
                        len = sqrt(sum(ab.^2));
                        lenkeeps = len~=0;
                        ab(:,lenkeeps) = ab(:,lenkeeps) ./ repmat(len(lenkeeps),3,1);
                        adjbvec = ab;
                        
                        adjbvec(1) = -adjbvec(1);

                        % Signal equation with adjected bvec and bval  
                        Lest_dwi(i,j,k,v) = b0_vol(i,j,k)*exp(-1*adjbval*adjbvec'*DT_mat(:,:)*adjbvec) ;
		            	NLest_dwi(i,j,k,v) = abs(Lest_dwi(i,j,k,v) + real_noise(i,j,k,v) + img_noise(i,j,k,v)) ;
                        % Compute ADC 
                        %ADC = (log(Lest_dwi(i,j,k,v) / (b0_vol(i,j,k))) * (1 / (bv_b0 - ob)));
                        %fprintf('ADC corpt %f for volume %i\n',[ADC, v]);
                    end
                end
            end
        end
    end  
    save('noise1.mat','real_noise','img_noise')
    % Check if intended SNR is the actual SNR
    diff = Nest_dwi - est_dwi;
    diff = diff(mask_vol);
    intended_snr =  mean_signal / std(diff)
    diff = NLest_dwi - Lest_dwi;
    diff = diff(mask_vol);
    intended_snr =  mean_signal / std(diff)

    % Save signal esimated. b0 volume is added since it was not used for the computation
    %dwmri_est = zeros(size(dwmri_vols));
    %dwmri_est(:,:,:,1) = b0_vol ;
    %dwmri_est(:,:,:,2:end) = est_dwi ; 
    %nii = load_untouch_nii(dwi_path);
    %nii.img = dwmri_est;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_est_sig']),nii,'double');
    
    % Difference in original and signal esimation
    %diff = abs(est_dwi - dwi_vols);
    %nii.img = diff;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_est_org_diff']),nii,'double');
    %noise = std(diff,0,4);
    %snr_est = nanmean(est_dwi,4) ./ noise;
    %nanmean(snr_est(mask_vol));
   
    % Save the corrput signal estimated. b0 volume is added since it was not used for the computation
    %dwmri_Lest = zeros(size(dwmri_vols));
    %dwmri_Lest(:,:,:,1) = b0_vol ;
    %dwmri_Lest(:,:,:,2:end) = Lest_dwi ;
    %nii.img = dwmri_Lest;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_sig']),nii,'double');
    %noise = std(Lest_dwi,0,4);
    %snr_Lest = nanmean(Lest_dwi,4) ./ noise;

    % Difference in orginal signal and after corrput signal esimation
    %Ldiff = abs(Lest_dwi - dwi_vols);
    %nii.img = Ldiff;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_org_diff']),nii,'double');
    %Lnoise = std(Ldiff,0,4);
    %snr = nanmean(est_dwi,4) ./ Lnoise;
    %nanmean(snr(mask_vol))
    
    % Difference in before and after corrput signal esimation
    %Ldiff_est = abs(est_dwi - Lest_dwi);
    %nii.img = Ldiff_est;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_est_diff']),nii,'double');
    %Lnoise = std(Ldiff,0,4);
    %snr = nanmean(est_dwi,4) ./ Lnoise;
    %nanmean(snr(mask_vol))

    dwmri_Nest = zeros(size(dwmri_vols));
    dwmri_Nest(:,:,:,1) = b0_vol ;
    dwmri_Nest(:,:,:,2:end) = Nest_dwi ;
    nii = load_untouch_nii(dwi_path);
    nii.img = dwmri_Nest;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_est_sig']),nii,'double');

    dwmri_NLest = zeros(size(dwmri_vols));
    dwmri_NLest(:,:,:,1) = b0_vol ;
    dwmri_NLest(:,:,:,2:end) = NLest_dwi ;
    nii = load_untouch_nii(dwi_path);
    nii.img = dwmri_NLest;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_sig']),nii,'double');
   
end
