function axial_d = compute_AD(D)
    % compute_PEV Compute PEV of diffusion tensor
    % Inputs
    % D         diffusion tensor
    % Output
    % PEV        PEV of diffusion tensor
    % v = eigen vector e = eigen values
    [v, e] = eig(D);
    e = diag(e);
    [max_eig, pos] = max(e);
    
    % max eigen value's corresponding eigen vector is the primary direction
    axial_d(:) = max_eig;
end
