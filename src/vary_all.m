%function sim_dt_FA_MD(MD,FA)
min = 0;
FA_max = 0.9;
n = 100;
FA = min + rand(1,n)*(FA_max-min);

    MD = 0.0008;
    %FA = 0.3;
    FA_sim_corpt_x = zeros(size(1:n));
    FA_true_x = zeros(size(1:n));
    FA_sim_x = zeros(size(1:n));
    FA_corr_bx_x = zeros(size(1:n));
    FA_corr_sm_x = zeros(size(1:n));
    %MD = 0.0008;
    % Trace fixed
    %TrD = 0.0021; %derek jones
    %MD = TrD/3;
    
    % FA change
    %FA = 0.5;
    %l1 = (MD) * (1 + (2*FA) / sqrt( 3 - (2 * FA^2)));
    %lprep = (MD) * (1 - FA / sqrt( 3 - (2 * FA^2)));
    for j = 1:n
    lprep = MD * (1 - FA(j)/ (3-2*FA(j)^2)^(1/2));
    % for PVE
    %for i = 1:n
     
    %a1 = 0.0;
    %b1 = 0.1;
    %sz = [1 n];
    min = 0;
    theta_max = 100;
    phi_max = 90;

    theta = min + rand(1,n)*(theta_max-min); %random 0-180 2 * pi * unifrnd(a1,b1,sz);
    phi = min + rand(1,n)*(phi_max-min); % 0-90 acos(1 - 2 * unifrnd(a1,b1,sz));
        
    for i = 1:n
        x = sin(phi(i)) * cos(theta(i));
        y = sin(phi(i)) * sin(theta(i));
        z = cos(phi(i));
        PVE = [x y z];
        
        DT_true = sim_dt(PVE, MD,  lprep);
	

        
        [S_corpt, S, abvec, abval, g, b] = sim_signal(DT_true, L_mat);
         
         FA_true = compute_FA(DT_true);
         FA_true_x(i) = FA_true;
         
        [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
        FA_sim_corpt = compute_FA(D_sim_corpt);
        FA_sim_corpt_x(i) =  FA_sim_corpt;
        
        [D_corr_bx, exitcode] =  linear_vox_fit(1,S_corpt,abvec, abval);
        FA_corr_bx = compute_FA(D_corr_bx);
        disp(FA_corr_bx)
        FA_corr_bx_x(i) = FA_corr_bx;

        new_dwi_signal = correct_signal_sm(S_corpt,b);
        [D_corr_sm, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
        FA_corr_sm = compute_FA(D_corr_sm);
        FA_corr_sm_x(i) =  FA_corr_sm;

        [D_sim, exitcode] = linear_vox_fit(1,S,g, b);
        FA_sim = compute_FA(D_sim);
        FA_sim_x(i) = FA_sim;        
        
        %FA = compute_FA(D);
        %plotDTI(D_true,0.002);
    end
end
