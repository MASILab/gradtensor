function [DT_vol, exitcode_vol] = linear_vol_fit(b0_vol,dwi_vol,bvecs,bvals,mask_vol)
    % Performs linear DTI fit on a volume. Note that only symmetric 
    % constraint is applied. Positive definiteness is not enforced.
    %
    % INPUTS:
    %   b0_vol is the non-diffusion weighted volume
    %   dwi_vol are the diffusion weighted volumes. Assumes diffusion
    %       weightings are done along the 4th dimension
    %   Assumes size(bvecs) = [3 size(dwi_vol,4)]. Also assumes bvecs are
    %       with respect to voxels.
    %   Assumes size(bvals) = [1 size(dwi_vol,4)]
    %   mask_vol is optional
    %   
    % OUTPUTS:
    %   DT_vol has format of [Dxx Dxy Dxz Dyy Dyz Dzz] along 4th dimension
    %   exitcode_vol:
    %   exitcode of -1 = background data
    %   exitcode of 0 = fit was successful
    %   exitcode of 1 = non-finite value in diffusion tensor

    if ~exist('mask_vol','var') || isempty(mask_vol)
        mask_vol = true(size(b0_vol));
    end

%     good_vals = squeeze(dwi_vol ~= 0);
    
%     if sum(good_vals) > 0
%         bvecs = bvecs(:, good_vals);
%         bvals = bvals(:, good_vals);
%         dwi_vol = dwi_vol(good_vals);
%     end
    bad_vals = squeeze(dwi_vol == 0);
    dwi_vol(bad_vals) = 0.0005;
    
    % Compute A: -b*g*g'
    gx = bvecs(1,:)';
    gy = bvecs(2,:)';
    gz = bvecs(3,:)';
    A = -bsxfun(@times,bvals',[gx.^2 2*gx.*gy 2*gx.*gz gy.^2 2*gy.*gz gz.^2]);

    % Reshape to make DT fit a vector equation
    dwi_size = size(dwi_vol);
    b0_vol = reshape(b0_vol,[],1);
    dwi_vol = reshape(dwi_vol,[],dwi_size(4));

    % Fit using linear least squares
    %disp(size(pinv(A)));
    %disp(size(permute(log(abs(bsxfun(@rdivide,dwi_vol,b0_vol))),[2 1])));
    DT_vol = pinv(A)*permute(log(abs(bsxfun(@rdivide,dwi_vol,b0_vol))),[2 1]);
    DT_vol = reshape(permute(DT_vol,[2 1]),dwi_size(1),dwi_size(2),dwi_size(3),[]);

    % Clear out values outside mask (don't worry about using the mask 
    % during computation; this fit is fast without it since it is 
    % completely vectorized).
    DT_vol(~repmat(mask_vol,1,1,1,6)) = 0;

    % Set exit code volume
    exitcode_vol = -ones(size(mask_vol));
    exitcode_vol(mask_vol) = 0;
    exitcode_vol(any(~isfinite(DT_vol),4)) = 1;
end
