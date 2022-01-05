function MD = compute_MD(D)
    % compute_MD Compute MD of diffusion tensor
    % Inputs
    % D         diffusion tensor
    % Output
    % MD        MD of diffusion tensor
    % v = eigen vector e = eigen values
    [v, e] = eig(D);
    e = diag(e);
    MD = (e(1) + e(2) + e(3))./3;
end
