% FA/MD/PEV after and before gradient non-linear fields in Tensor simulation 
MD = 0.0008;
n = 1;

% Declare array
FA_sim_corpt_x = zeros(n*n*n*n,1);
FA_true_x = zeros(n*n*n*n,1);
FA_sim_x = zeros(n*n*n*n,1);
FA_corr_bx_x = zeros(n*n*n*n,1);
FA_corr_sm_x = zeros(n*n*n*n,1);
FA_corr_sm_fy_x = zeros(n*n*n*n,1);

MD_sim_corpt_x = zeros(n*n*n*n,1);
MD_true_x = zeros(n*n*n*n,1);
MD_sim_x = zeros(n*n*n*n,1);
MD_corr_bx_x = zeros(n*n*n*n,1);
MD_corr_sm_x = zeros(n*n*n*n,1);
MD_corr_sm_fy_x = zeros(n*n*n*n,1);

PEV_sim_corpt_x = zeros(n*n*n*n,3);
PEV_true_x = zeros(n*n*n*n,3);
PEV_sim_x = zeros(n*n*n*n,3);
PEV_corr_bx_x = zeros(n*n*n*n,3);
PEV_corr_sm_x = zeros(n*n*n*n,3);
PEV_corr_sm_fy_x = zeros(n*n*n*n,3);

FA_all = zeros(1,n*n*n*n);
phi_all = zeros(n*n*n*n,1);
theta_all = zeros(n*n*n*n,1);
Ldet_all = zeros(n*n*n*n,1);

% Generate evenly distributed FA
FA = linspace(0,1,n);
%lprep = zeros(1,n);
for i = 1:n
   lprep(i) = MD * (1 - FA(i)/ (3-2*FA(i)^2)^(1/2));
end

% Rank-ordered expected LR field and selected 100 by choosing the first through 100th percentiles
[L_det, vL] = vary_L_det();
A = L_det;
blah = A(A>0);
[out,idx] = sort(blah,'ascend');
indices_to_chose = round(linspace(1,626688,n));
vec = out(indices_to_chose);
%det_vals = zeros(1,n);
locations = zeros(3,n);
for k = 1:n
  [r,c,v] = ind2sub(size(A),find(A == vec(k)));
  det_vals(k) = A(r,c,v);
  locations(:,k) = [r c v];
end

% Theta and phi
theta_ls = linspace(0, 2*pi, n);
phi_ls = linspace(0, acos(1 - 2*0.5), n);
PVE = zeros(n*n,3);
idx = 0;

% Compute for every theta, phi, LR and FA 
for i = 1:n
    for j = 1:n
        theta = theta_ls(i);
        phi = phi_ls(j);
        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);
        PVE(( (i-1) * n) + j,:) = [x y z];
        for p = 1:n
            for o = 1:n
                % Saving all the variables when looped
                idx = idx + 1;
                theta_all(idx) =  theta;
                phi_all(idx) = phi;
                FA_all(idx) =  FA(p);
                Ldet_all(idx) =  det_vals(o);

                % Get L matrix at position xyz
                xyz = locations(:,o);
                L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));

                % Simulate the tensor for given PEV, FA
                DT_true = ten_compute_dt([x y z], MD,  lprep(p));

                % Simulate signal for fake tensor with L mat 
                [S_corpt, S, abvec, abval, g, b] = ten_compute_corrput_signal(DT_true, L_mat);

                % Compute FA + MD + PEV before corrpution 
                FA_true = compute_FA(DT_true);
                 FA_true_x(idx) =  FA_true;
                 MD_true = compute_MD(DT_true);
                 MD_true_x(idx) = MD_true;
                 PEV_true = compute_primary_eigvec(DT_true);
                 PEV_true_x(idx,:) =  PEV_true;

                % DT fit corpt signal
                [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
                % Compute FA + MD + PEV after corrpution
		FA_sim_corpt = compute_FA(D_sim_corpt);
                 FA_sim_corpt_x(idx) =   FA_sim_corpt;
                 MD_sim_corpt = compute_MD(D_sim_corpt);
                 MD_sim_corpt_x(idx) =   MD_sim_corpt;
                 PEV_sim_corpt = compute_primary_eigvec(D_sim_corpt);
                 PEV_sim_corpt_x(idx,:) =   PEV_sim_corpt;

                %plotDTI(DT_true, 0.002, [0.25, 0.33, 0.40]);
            end
        end
    end
end


