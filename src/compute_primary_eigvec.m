function primary_vec_vol = compute_PEV(D)
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
    primary = v(:,pos);
    eig_vol(:) = e;
    primary_vec_vol(:) = primary;
end
