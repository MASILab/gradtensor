function [DT_mat, exitcode] = linear_vox_fit(b0,dwi,bvecs,bvals)
    % Performs linear DTI fit on a voxel. Note that only symmetric 
    % constraint is applied. Positive definiteness is not enforced.
    %
    % INPUTS:
    %   b0 is the non-diffusion weighted value
    %   dwi are diffusion weighted values
    %   Assumes size(bvecs) = [3 length(dwi)]. Also assumes bvecs are
    %       with respect to voxels.
    %   Assumes size(bvals) = [1 length(dwi)]
    %   
    % OUTPUTS:
    %   DT_mat has format of [Dxx Dxy Dxz 
    %                         Dxy Dyy Dyz 
    %                         Dxz Dyz Dzz]
    %   exitcode info in vol_fit

    % Just reshape and call the vol_fit since it is already vectorized.
    % Note that exitcode is set in vol_fit.
    [DT, exitcode] = linear_vol_fit(b0,reshape(dwi,1,1,1,length(dwi)),bvecs,bvals);

    % Form DT_mat from DT
    DT_mat = zeros(3);
    DT_mat(1,:) = DT(1:3);
    DT_mat(2,2:end) = DT(4:5);
    DT_mat(3,3) = DT(6);
    % Copy diagonal elements
    DT_mat(2,1) = DT_mat(1,2);
    DT_mat(3,1:2) = DT_mat(1:2,3);
end
