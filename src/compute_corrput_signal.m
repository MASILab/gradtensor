function compute_corrput_signal(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name, L_path, org_bvec_path, org_bval_path)
    % compute_corrput_signal Computes signal with LR induced in it 
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
    % Load the original bvec and bval 
    org_bvecs = importdata(org_bvec_path);
    org_bvals = importdata(org_bval_path);
    
    % dwi and non-dwi
    ind_b0 = find(~org_bvals);
    ind_non_b0 = find(org_bvals);
    
    org_bvecs = org_bvecs(:,ind_non_b0);
    org_bvals = org_bvals(ind_non_b0);
    
    all_b0_vol = dwmri_vols(:,:,:,ind_b0);
    b0_vol = mean(all_b0_vol,4);
    dwi_vols = dwmri_vols(:,:,:,ind_non_b0);

    % Load resampled L
    % Resample grad tensor to DTI image space
    % Unzip
    if strcmp(L_path(end-2:end),'.gz')
        system(['gunzip -kf ' L_path]);
        L_path = L_path(1:end-3);
    end
    flags = struct( ...
        'mask',true, ...
        'mean',false, ...
        'interp',1, ...
        'which',1, ...
        'wrap',[0 0 0], ...
        'prefix','r' ...
        );
    spm_reslice({dwi_path; L_path},flags);
    [p,n,e] = fileparts(L_path);
    rLimg_file = fullfile(p,['r' n e]);
    movefile(rLimg_file,fullfile(out_dir,'L_resamp.nii'));
    rLimg_file = fullfile(out_dir,'L_resamp.nii');

    VL = spm_vol(rLimg_file);
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
        disp('making mask')
        mask_vol = true(size(b0_vol));
    else
        %gunzip(mask_path)
        mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
        mask_vol = logical(mask_vol);
    end
    
    % Volume 2 to end for bvec and bval
    bvec_vols = bvec_vols(:,:,:,2:end,:);
    bval_vols = bval_vols(:,:,:,2:end);

    % Initialize DT and exitcode volumes
    exitcode_vol = zeros(size(b0_vol));
    eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
    primary_vec_vol = zeros(size(eig_vol)); 
    est_dwi = zeros(size(dwi_vols));
    Lest_dwi = zeros(size(dwi_vols));

    % Cycle over and compute DT voxel-wise
    for i = 1:size(dwi_vols,1)
        for j = 1:size(dwi_vols,2)
            for k = 1:size(dwi_vols,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals               
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';
                    %bvecs = squeeze(bvec_vols(i,j,k,:,:))';
                    %bvals = squeeze(bval_vols(i,j,k,:))';
                    L_mat = squeeze(vL(:,:,i,j,k));

                    % Get linear model - use the corrected bvec and bval - considered as ground 
                    % truth tensor
                    %[DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvecs,bvals);
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,org_bvecs,org_bvals);
                   
                    % Induce corrpution (script adpated from compute_b_images.m) and signal equation
                    for v = 1:length(org_bvals)
                        %g = bvecs(:,v);
                        %b = bvals(v);
                        
                        og = org_bvecs(:,v);
                        ob = org_bvals(v);
                        
                        % DWI signal with no corput
                        est_dwi(i,j,k,v) =  b0_vol(i,j,k)*exp(-1*ob*og'*DT_mat(:,:)*og);

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
                        Lest_dwi(i,j,k,v) = b0_vol(i,j,k)*exp(-1*adjbval*adjbvec'*DT_mat(:,:)*adjbvec);

                        % Compute ADC 
                        %ADC = (log(Lest_dwi(i,j,k,v) / (b0_vol(i,j,k))) * (1 / (bv_b0 - ob)));
                        %fprintf('ADC corpt %f for volume %i\n',[ADC, v]);
                    end
                end
            end
        end
    end  

    % Save signal esimated. b0 volume is added since it was not used for the computation
    dwmri_est = zeros(size(dwmri_vols));
    dwmri_est(:,:,:,ind_b0) = all_b0_vol ;
    dwmri_est(:,:,:,ind_non_b0) = est_dwi ; 
    nii = load_untouch_nii(dwi_path);
    nii.img = dwmri_est;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_est_sig']),nii,'double');
    
    % Difference in original and signal esimation
    %diff = abs(est_dwi - dwi_vols);
    %nii.img = diff;
    %nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_est_org_diff']),nii,'double');
    %noise = std(diff,0,4);
    %snr = nanmean(est_dwi,4) ./ noise;
    %nanmean(snr(mask_vol))
   
    % Save the corrput signal estimated. b0 volume is added since it was not used for the computation
    dwmri_Lest = zeros(size(dwmri_vols));
    dwmri_Lest(:,:,:,ind_b0) = all_b0_vol ;
    dwmri_Lest(:,:,:,ind_non_b0) = Lest_dwi ;
    nii.img = dwmri_Lest;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lest_sig']),nii,'double');

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
    
end
