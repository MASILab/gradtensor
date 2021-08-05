function primary_vec_vol = compute_PEV(D)
    [v, e] = eig(D);
    e = diag(e);
    [max_eig, pos] = max(e);
    primary = v(:,pos);
    eig_vol(:) = e;
    primary_vec_vol(:) = primary;
end
