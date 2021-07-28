function FA = compute_FA(D)
    [v, e] = eig(D);
    e = diag(e);
    FA = sqrt(1/2) .* (sqrt( (e(1) - e(2)).^2 + (e(2) - e(3)).^2  + (e(3) - e(1)).^2 ) ./ sqrt(e(1).^2 + e(2).^2 + e(3).^2));
end
