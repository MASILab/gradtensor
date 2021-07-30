function MD = compute_MD(D)
    [v, e] = eig(D);
    e = diag(e);
    MD = (e(1) + e(2) + e(3))./3;
end
