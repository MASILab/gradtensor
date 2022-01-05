function FA = compute_FA(D)
    % compute_FA Compute FA of diffusion tensor
    % Inputs 
    % D 	diffusion tensor
    % Output
    % FA 	FA of diffusion tensor
    % v = eigen vector e = eigen values
    [v, e] = eig(D);
    e = diag(e);
    FA = sqrt(1/2) .* (sqrt( (e(1) - e(2)).^2 + (e(2) - e(3)).^2  + (e(3) - e(1)).^2 ) ./ sqrt(e(1).^2 + e(2).^2 + e(3).^2));
end
