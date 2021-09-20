MD = 0.0008;
n = 40;

% declare array
FA_sim_corpt_x = [];
FA_true_x = [];
FA_sim_x = [];
FA_corr_bx_x = [];
FA_corr_sm_x = [];
FA_corr_sm_fy_x =[];

MD_sim_corpt_x = [];
MD_true_x = [];
MD_sim_x = [];
MD_corr_bx_x =[];
MD_corr_sm_x = [];
MD_corr_sm_fy_x = [];

PEV_sim_corpt_x = [];
PEV_true_x = [];
PEV_sim_x = [];
PEV_corr_bx_x = [];
PEV_corr_sm_x = [];
PEV_corr_sm_fy_x = [];

FA_all = [];
phi_all = [];
theta_all = [];
Ldet_all = [];

% generate evenly distributed FA
FA = linspace(0,1,n);
%lprep = zeros(1,n);
for i = 1:n
   lprep(i) = MD * (1 - FA(i)/ (3-2*FA(i)^2)^(1/2));
end

% rank-ordered expected LR field and selected 100 by choosing the first through 100th percentiles
[L_det, vL] = vary_L_det();
A = L_det;
blah = A(A>0);
[out,idx] = sort(blah,'ascend');
indices_to_chose = round(linspace(1,626688,n));
vec = out(indices_to_chose);
%det_vals = zeros(1,n);
%locations = zeros(3,n);
for k = 1:n
  [r,c,v] = ind2sub(size(A),find(A == vec(k)));
  det_vals(k) = A(r,c,v);
  locations(:,k) = [r c v];
end

% theta and phi
theta_ls = linspace(0, 2*pi, n);
phi_ls = linspace(0, acos(1 - 2*0.5), n);
%PVE = zeros(n*n,3);
lala = 0;

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
                 % saving all the variables when looped
                 % n((i-1) + (j-1) + (p-1)) + o
                %idx = n * (i + j + p - 3) + o;
                theta_all = [theta_all; theta];
                phi_all = [phi_all; phi];
                FA_all =  [FA_all; FA(p)];
                Ldet_all = [Ldet_all; det_vals(o)];

                % get L matrix at position xyz
                xyz = locations(:,o);
                %disp(xyz);
                L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));
                %L_mat(:,:) = vL(:,:,58,52,41);

                % simulate the tensor for given PEV, FA
                DT_true = sim_dt([x y z], MD,  lprep(p));

                % simulate signal for fake tensor with L mat 
                [S_corpt, S, abvec, abval, g, b] = sim_signal(DT_true, L_mat);

                % compute FA + MD + PEV before corrpution 
                FA_true = compute_FA(DT_true);
                %disp(FA_true)
                FA_true_x =  [FA_true_x; FA_true];
                MD_true = compute_MD(DT_true);
                MD_true_x = [MD_true_x; MD_true];
                PEV_true = compute_primary_eigvec(DT_true);
                PEV_true_x =  [PEV_true_x; PEV_true];

                % DT fit corpt signal
                [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
                % compute FA + MD + PEV after corrpution
                FA_sim_corpt = compute_FA(D_sim_corpt);
                FA_sim_corpt_x = [FA_sim_corpt_x;  FA_sim_corpt];
                MD_sim_corpt = compute_MD(D_sim_corpt);
                MD_sim_corpt_x =  [MD_sim_corpt_x;  MD_sim_corpt];
                PEV_sim_corpt = compute_primary_eigvec(D_sim_corpt);
                PEV_sim_corpt_x = [PEV_sim_corpt_x;  PEV_sim_corpt];

                % DT fit baxter corrected with adjected bvec and adjected bval
                [D_corr_bx, exitcode] =  linear_vox_fit(1,S_corpt,abvec, abval);
                % FA + MD + PVE baxter correction
                FA_corr_bx = compute_FA(D_corr_bx);
                FA_corr_bx_x =  [FA_corr_bx_x; FA_corr_bx];
                MD_corr_bx = compute_MD(D_corr_bx);
                MD_corr_bx_x =  [MD_corr_bx_x; MD_corr_bx];
                PEV_corr_bx = compute_primary_eigvec(D_corr_bx);
                PEV_corr_bx_x = [PEV_corr_bx_x; PEV_corr_bx];

                % compute new signal with det + DT fit 
%                 new_dwi_signal = correct_signal_sm(S_corpt,g,b,L_mat);
%                 [D_corr_sm, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
%                 % FA + MD + PVE diffusivity method
%                 FA_corr_sm = compute_FA(D_corr_sm);
%                 FA_corr_sm_x =  [FA_corr_sm_x; FA_corr_sm];
%                 MD_corr_sm = compute_MD(D_corr_sm);
%                 MD_corr_sm_x =  [MD_corr_sm_x; MD_corr_sm];
%                 PEV_corr_sm = compute_primary_eigvec(D_corr_sm);
%                 PEV_corr_sm_x =  [PEV_corr_sm_x; PEV_corr_sm];

                %compute new signal with scaling the bvec + DT fit
                fy_new_dwi_signal = bl_correct_signal_sm(S_corpt,g,b,abvec,abval,L_mat);
                [D_corr_sm_fy, exitcode] = linear_vox_fit(1,fy_new_dwi_signal,g, b);
                % FA + MD + PVE anisotropy method
                FA_corr_sm_fy = compute_FA(D_corr_sm_fy);
                FA_corr_sm_fy_x = [FA_corr_sm_fy_x;  FA_corr_sm_fy];
                MD_corr_sm_fy = compute_MD(D_corr_sm_fy);
                MD_corr_sm_fy_x = [MD_corr_sm_fy_x;  MD_corr_sm_fy];
                PEV_corr_sm_fy = compute_primary_eigvec(D_corr_sm_fy);
                PEV_corr_sm_fy_x =  [PEV_corr_sm_fy_x ; PEV_corr_sm_fy];

                % DT for simulate tensor without corpt % should be like my true
                % tensor
%                 [D_sim, exitcode] = linear_vox_fit(1,S,g, b);
%                 FA_sim = compute_FA(D_sim);
%                 FA_sim_x = [FA_sim_x; FA_sim];

                %plotDTI(DT_true, 0.002)%,'b');
                lala = lala + 1;
                %disp(lala);
            end
        end
    end
end


