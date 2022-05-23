function radial_d = compute_RD(D)
    % compute_PEV Compute PEV of diffusion tensor
    % Inputs
    % D         diffusion tensor
    % Output
    % PEV        PEV of diffusion tensor
    % v = eigen vector e = eigen values
    [v, e] = eig(D);
    e = diag(e);
    [max_eig, pos] = max(e);

    % sum of the least two 
    [min_eig, ter_pos] = min(e);
    possible_positions = [1 2 3];
    not_max = possible_positions(possible_positions~=pos);
    sec_pos = not_max(not_max~=ter_pos);
    sec_eig = e(sec_pos);
    radial_d(:) = (min_eig + sec_eig) / 2;
end
