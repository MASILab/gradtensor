function voxel_to_sig(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name,bval_file,bvec_file)
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
    dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
    bvec = load(bvec_file);
    bval = load(bval_file);
    bval = bval';

    % dwi and non-dwi
    ind_b0 = find(~bval);
    ind_non_b0 = find(bval);

    all_b0_vol = dwmri_vols(:,:,:,ind_b0);
    b0_vol = mean(all_b0_vol,4);
    dwi_vols = dwmri_vols(:,:,:,ind_non_b0);
    size(dwi_vols)
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
    corr_dwi = zeros(size(dwi_vols));
    size(corr_dwi)
    % Cycle over and compute DT voxel-wise
    for i = 1:size(dwi_vols,1)
        for j = 1:size(dwi_vols,2)
            for k = 1:size(dwi_vols,3)
                if mask_vol(i,j,k)
                    % Get b0, dwi, bvecs and bvals
                    b0 = b0_vol(i,j,k);
                    dwi = squeeze(dwi_vols(i,j,k,:))';
                    bvecs = squeeze(bvec_vols(i,j,k,:,:))';
                    bvals = squeeze(bval_vols(i,j,k,:))';

                    % Get linear model
                    [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvecs,bvals); 
                    for v = 1:size(dwi_vols,4)
                        %g = bvecs(:,v);
                        %b = bvals(v);

                        og = bvec(:,v);
                        ob = bval(v);

                        % Est signal from voxelwise DT mat
                        corr_dwi(i,j,k,v) =  b0_vol(i,j,k)*exp(-1*ob*og'*DT_mat(:,:)*og);
                     
                  
                    end
                end
            end
        end
    end

    % Compute and save the MD, FA, PEV of the DTI
    dwmri_corr = zeros(size(dwmri_vols));
    dwmri_corr(:,:,:,ind_b0) = all_b0_vol ;
    dwmri_corr(:,:,:,ind_non_b0) = corr_dwi ;
    nii = load_untouch_nii(dwi_path);
    dwmri_corr(isnan(dwmri_corr)) = 0;
    nii.img = dwmri_corr;
    nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_Lcorrected_sig']),nii,'double');
    
end
