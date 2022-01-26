% FA/MD/PEV after and before gradient non-linear fields in Tensor simulation
MD = 0.0008;
n = 100;
nv = 19406;

% Declare array
FA_sim_corptn_x = zeros(nv*n*n,1);
FA_sim_corpt_x = zeros(nv*n*n,1);
FA_true_x = zeros(nv*n*n,1);
FA_corr_bx_x = zeros(nv*n*n,1);
FA_corr_sm_fy_x = zeros(nv*n*n,1);

MD_sim_corptn_x = zeros(nv*n*n,1);
MD_sim_corpt_x = zeros(nv*n*n,1);
MD_true_x = zeros(nv*n*n,1);
MD_corr_bx_x = zeros(nv*n*n,1);
MD_corr_sm_fy_x = zeros(nv*n*n,1);

PEV_sim_corptn_x = zeros(nv*n*n,3);
PEV_sim_corpt_x = zeros(nv*n*n,3);
PEV_true_x = zeros(nv*n*n,3);
PEV_corr_bx_x = zeros(nv*n*n,3);
PEV_corr_sm_fy_x = zeros(nv*n*n,3);

FA_all = zeros(1,nv*n*n);
Ldet_all = zeros(nv*n*n,1);

% Generate evenly distributed FA
FA = linspace(0.01,1,n);
for i = 1:n
   lprep(i) = MD * (1 - FA(i)/ (3-2*FA(i)^2)^(1/2));
end

% Rank-ordered expected LR field and selected 100 by choosing the first through 100th percentiles
[L_det, vL] = vary_L_det();
A = L_det;
blah = A(A>0);
[out,index] = sort(blah,'ascend');
indices_to_chose = round(linspace(1,626688,n));
vec = out(indices_to_chose);
locations = zeros(3,n);
for k = 1:n
  [r,c,v] = ind2sub(size(A),find(A == vec(k)));
  det_vals(k) = A(r,c,v);
  locations(:,k) = [r c v];
end

% Theta and phi for dipy uniformly distributed sphere
load('theta_sphere100.mat');
load('phi_sphere100.mat');
load('vertices_sphere100.mat');
idx = 0;
for i = 1:nv
        particle_cord = vertices(i, :);
        for p = 1:n
            for o = 1:n
                 % saving all the variables when looped
                idx = idx + 1;
                FA_all(idx) =  FA(p);
                Ldet_all(idx) =  det_vals(o);

                % Get L matrix at position xyz
                xyz = locations(:,o);
                L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));

                % Simulate the tensor for given PEV, FA
                x = particle_cord(1);
                y = particle_cord(2);
                z = particle_cord(3);
                DT_true = ten_compute_dt([x y z], MD,  lprep(p));

                % Simulate signal for fake tensor with L mat 
                [S_corpt, S, abvec, abval, g, b] = ten_compute_corrput_signal100(DT_true, L_mat);

                % Compute FA + MD + PEV before corrpution 
                FA_true = compute_FA(DT_true);
                FA_true_x(idx) =  FA_true;
                MD_true = compute_MD(DT_true);
                MD_true_x(idx) = MD_true;
                PEV_true = compute_primary_eigvec(DT_true);
                PEV_true_x(idx,:) =  PEV_true;
  		
		% DT fit corpt signal w only noise
                [D_sim_corptn, exitcode] = linear_vox_fit(1,S,g, b);
                % Compute FA + MD + PEV after corrpution
                FA_sim_corptn = compute_FA(D_sim_corptn);
                FA_sim_corptn_x(idx) =   FA_sim_corptn;
                MD_sim_corptn = compute_MD(D_sim_corptn);
                MD_sim_corptn_x(idx) =   MD_sim_corptn;
                PEV_sim_corptn = compute_primary_eigvec(D_sim_corptn);
                PEV_sim_corptn_x(idx,:) =   PEV_sim_corptn;

                % DT fit corpt signal
                [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
                % Compute FA + MD + PEV after corrpution
                FA_sim_corpt = compute_FA(D_sim_corpt);
                FA_sim_corpt_x(idx) =   FA_sim_corpt;
                MD_sim_corpt = compute_MD(D_sim_corpt);
                MD_sim_corpt_x(idx) =   MD_sim_corpt;
                PEV_sim_corpt = compute_primary_eigvec(D_sim_corpt);
                PEV_sim_corpt_x(idx,:) =   PEV_sim_corpt;
		
                % DT fit full corrected with adjected bvec and adjected bval
                 [D_corr_bx, exitcode] =  linear_vox_fit(1,S_corpt,abvec, abval);
                 % compute_ADC_(D_corr_bx, g, b);
                 % FA + MD + PVE full correction method
                 FA_corr_bx = compute_FA(D_corr_bx);
                 FA_corr_bx_x(idx) =  FA_corr_bx;
                 MD_corr_bx = compute_MD(D_corr_bx);
                 MD_corr_bx_x(idx) =  MD_corr_bx;
                 PEV_corr_bx = compute_primary_eigvec(D_corr_bx);
                 PEV_corr_bx_x(idx,:) = PEV_corr_bx;

                 % Compute new signal with scaling the bvec + DT fit
                 new_dwi_signal = ten_compute_approx_corr(S_corpt,g,b,L_mat);
                 [D_corr_sm_fy, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
                 % FA + MD + PVE approximate correction method
                 FA_corr_sm_fy = compute_FA(D_corr_sm_fy);
                 FA_corr_sm_fy_x(idx) =   FA_corr_sm_fy;
                 MD_corr_sm_fy = compute_MD(D_corr_sm_fy);
                 MD_corr_sm_fy_x(idx) =   MD_corr_sm_fy;
                 PEV_corr_sm_fy = compute_primary_eigvec(D_corr_sm_fy);
                 PEV_corr_sm_fy_x(idx,:) =  PEV_corr_sm_fy;

%                 disp('--------FA-----------')
%                 disp(FA_true_x)
% 		         disp(FA_sim_corpt_x)
% 		         disp(FA_corr_bx_x)
% 		         disp(FA_corr_sm_fy_x)
%                 disp('---------MD-----------')
%                 disp(MD_true_x)
% 		         disp(MD_sim_corpt_x)
% 		         disp(MD_corr_bx_x)
% 		         disp(MD_corr_sm_fy_x)
%                 disp('---------PEV----------')
%                 disp(PEV_true_x)
% 		         disp(PEV_sim_corpt_x)
% 		         disp(PEV_corr_bx_x)
% 		         disp(PEV_corr_sm_fy_x)
                %plotDTI(D_corr_bx, 0.002, [0.4660, 0.6740, 0.1880]);%[0.4940, 0.1840, 0.5560]);
            end
        end
end
save('/nfs/masi/kanakap/CorrectionLReffectsnoise100.mat','-v7.3');



