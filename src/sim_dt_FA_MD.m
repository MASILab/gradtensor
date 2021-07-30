function sim_dt_FA_MD(MD,FA)
    %MD = 0.0008;
    % Trace fixed
    %TrD = 0.0021; %derek jones
    %MD = TrD/3;

    % FA change
    %FA = 0.5;
    %l1 = (MD) * (1 + (2*FA) / sqrt( 3 - (2 * FA^2)));
    %lprep = (MD) * (1 - FA / sqrt( 3 - (2 * FA^2)));
    lprep = MD * (1 - FA/ (3-2*FA^2)^(1/2));

% for PVE
    n = 100;
    for i = 1:n
        PVE = [0, 0, 0];

        while norm(PVE) < 0.001
            x = randn();
            y = randn();
            z = randn();
            PVE = [x y z];
        end
        PVE = PVE / norm(PVE);
        disp(PVE)
        D = sim_dt(PVE, MD,  lprep);
        FA = compute_FA(D);
        plotDTI(D,0.002);
    end
end