function dti_voxel_fit_masiver(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name,bval_file,bvec_file)
    % dti_voxel_fit Performs linear DTI fit given input bfields. Note that only symmetric 
    % constraint is applied. Positive definiteness is not enforced. 
    %
    % Inputs
    %   dwi_path 	Path to the diffusion weighted nifti. Assumes b0 at first volume
    %   bvec_paths 	Cell array of paths to bvec niftis. Assumes these are with respect to voxels.
    %   bval_paths 	Ccell array of paths to bval niftis.
    %   mask_path 	Optional path to mask nifti.
    %   out_dir 	Path to a directory in which to save metrics
    %   out_name 	Prefix to the generated metric nifti filenames
    
    % Load data
    dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
    bvec = load(bvec_file);
    bval = load(bval_file);

    % dwi and non-dwi
    ind_b0 = find(~bval);
    ind_non_b0 = find(bval);

    all_b0_vol = dwi_vols(:,:,:,ind_b0);
    b0_vol = mean(all_b0_vol,4);
    dwi_vols = dwi_vols(:,:,:,ind_non_b0);
    
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
    if ~exist('mask_path','var') || isempty(mask_path)
        mask_vol = true(size(b0_vol));
    else
        %gunzip(mask_path)
        mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
        mask_vol = logical(mask_vol);
    end
    
    bvec_vols = bvec_vols(:,:,:,ind_non_b0,:);
    bval_vols = bval_vols(:,:,:,ind_non_b0);

    % Initialize DT and exitcode volumes
    exitcode_vol = zeros(size(b0_vol));
    eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
    primary_vec_vol = zeros(size(eig_vol)); 
    AD = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3));
    RD = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3));

    % Cycle over and compute DT voxel-wise
    for i = 1:size(mask_vol,1)
        for j = 1:size(mask_vol,2)
            for k = 1:size(mask_vol,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';
                    bvecs = squeeze(bvec_vols(i,j,k,:,:))';
                    bvals = squeeze(bval_vols(i,j,k,:))';

                    % Get linear model
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvecs,bvals); 
                     
                    if sum(isnan(DT_mat(:))) == 0 && sum(isinf(DT_mat(:))) == 0
                        %DT_rot_lest_corr_bx_posB(:,:,i,j,k) = DT_mat(:,:);
                    	[v, e] = eig(DT_mat);
                    	e = diag(e);
                        if max(e) == min(e)
                            [max_eig, pos] = max(e);
                            
                            min_eig = e(2);
                            sec_eig = e(3);
                            
                            primary = v(:,pos);	
                            eig_vol(i,j,k,:) = e;
                            primary_vec_vol(i,j,k,:) = primary;

                            AD(i,j,k) = max_eig;
                            RD(i,j,k) = (min_eig + sec_eig) / 2;
                        else
                            [max_eig, pos] = max(e);
                            [min_eig, ter_pos] = min(e);
                            possible_positions = [1 2 3];
                            not_max = possible_positions(possible_positions~=pos);
                            sec_pos = not_max(not_max~=ter_pos);
                            sec_eig = e(sec_pos);

                            primary = v(:,pos);	
                            eig_vol(i,j,k,:) = e;
                            primary_vec_vol(i,j,k,:) = primary;

                            AD(i,j,k) = max_eig;
                            RD(i,j,k) = (min_eig + sec_eig) / 2;
                        end
                        exitcode_vol(i,j,k) = exitcode;
                    end
                end
            end
        end
    end

    % Compute and save the MD, FA, PEV of the DTI
    MD = (eig_vol(:,:,:,1) + eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./3;
    FA = sqrt(1/2) .* (sqrt( (eig_vol(:,:,:,1) - eig_vol(:,:,:,2)).^2 + (eig_vol(:,:,:,2) - eig_vol(:,:,:,3)).^2  + (eig_vol(:,:,:,3) - eig_vol(:,:,:,1)).^2 ) ./ sqrt(eig_vol(:,:,:,1).^2 + eig_vol(:,:,:,2).^2 + eig_vol(:,:,:,3).^2));
    %AD = eig_vol(:,:,:,1);
    %RD = (eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./2;
    nii = load_untouch_nii(dwi_path);
    nii.img = MD;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_md']),nii,'double');
    nii.img = FA;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_fa']),nii,'double');
    nii.img = primary_vec_vol;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_primary_eigvec']),nii,'double');
       nii.img = AD;
     nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_ad']),nii,'double');
    nii.img = RD;
     nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_rd']),nii,'double'); 
    
end
