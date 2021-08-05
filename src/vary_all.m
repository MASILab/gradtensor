MD = 0.0008;
n = 100;

% declare array
FA_sim_corpt_x = [];
FA_true_x = [];
FA_sim_x = [];
FA_corr_bx_x = [];
FA_corr_sm_x = [];
FA_corr_sm_fy_x = [];

MD_sim_corpt_x = [];
MD_true_x = [];
MD_sim_x = [];
MD_corr_bx_x = [];
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
min_value = 0;
FA_max = 0.9;
FA = linspace(0,1,n);
%FA = min + rand(1,n)*(FA_max-min);
for i = 1:n
   lprep(i) = MD * (1 - FA(i)/ (3-2*FA(i)^2)^(1/2));
end

% evenly distributed theta and phi
theta_max = 180;
phi_max = 90;
theta = linspace(0,180,n);
phi = linspace(0,90,n);
%theta = min + rand(1,n)*(theta_max-min); %random 0-180 2 * pi * unifrnd(a1,b1,sz);
%phi = min + rand(1,n)*(phi_max-min); % 0-90 acos(1 - 2 * unifrnd(a1,b1,sz));
for j = 1:n
    x = sin(phi(j)) * cos(theta(j));
    y = sin(phi(j)) * sin(theta(j));
    z = cos(phi(j));
    PVE(j,:) = [x y z];
end

% rank-ordered expected LR field and selected 100 by choosing the first through 100th percentiles
[L_det, vL] = vary_L_det();
A = L_det;
blah = A(A>0);
[out,idx] = sort(blah,'ascend');
indices_to_chose = round(linspace(1,119862,n));
vec = out(indices_to_chose);
for k = 1:n
  [r,c,v] = ind2sub(size(A),find(A == vec(k)));
  det_vals(k) = A(r,c,v);
  locations(:,k) = [r c v];
end
lala = 0;
for i = 1:n
    for j = 1:n
        for o = 1:n
	    % saving all the variables when looped
            theta_all = [theta_all; theta(i)];
            phi_all = [phi_all; phi(i)];
            FA_all = [FA_all; FA(j)];
            Ldet_all = [Ldet_all; det_vals(o)];
            
	    % get L matrix at position xyz
            xyz = locations(:,o);
            disp(xyz);
            L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));
            
	    % simulate the tensor for given PEV, FA
            DT_true = sim_dt(PVE(i,:), MD,  lprep(j));
            
	    % simulate signal for fake tensor with L mat 
            [S_corpt, S, abvec, abval, g, b] = sim_signal(DT_true, L_mat);
		
	    % compute FA + MD + PEV before corrpution 
            FA_true = compute_FA(DT_true);
            FA_true_x = [FA_true_x; FA_true];
	    MD_true = compute_MD(DT_true);
            MD_true_x = [MD_true_x; MD_true];
	    PEV_true = compute_primary_eigvec(DT_true);
            PEV_true_x = [PEV_true_x; PEV_true];
	    
	    % DT fit corpt signal
            [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
	    % compute FA + MD + PEV after corrpution
            FA_sim_corpt = compute_FA(D_sim_corpt);
            FA_sim_corpt_x =  [FA_sim_corpt_x; FA_sim_corpt];
	    MD_sim_corpt = compute_MD(D_sim_corpt);
            MD_sim_corpt_x =  [MD_sim_corpt_x; MD_sim_corpt];
	    PEV_sim_corpt = compute_primary_eigvec(D_sim_corpt);
            PEV_sim_corpt_x =  [PEV_sim_corpt_x; PEV_sim_corpt];

	    % DT fit baxter corrected with adjected bvec and adjected bval
            [D_corr_bx, exitcode] =  linear_vox_fit(1,S_corpt,abvec, abval);
	    % FA + MD + PVE baxter correction
            FA_corr_bx = compute_FA(D_corr_bx);
            FA_corr_bx_x = [FA_corr_bx_x; FA_corr_bx];
	    MD_corr_bx = compute_MD(D_corr_bx);
            MD_corr_bx_x = [MD_corr_bx_x; MD_corr_bx];
	    PEV_corr_bx = compute_primary_eigvec(D_corr_bx);
            PEV_corr_bx_x = [PEV_corr_bx_x; PEV_corr_bx];
            
	    % compute new signal with det + DT fit 
            new_dwi_signal = correct_signal_sm(S_corpt,g,b,L_mat);
            [D_corr_sm, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
	    % FA + MD + PVE diffusivity method
            FA_corr_sm = compute_FA(D_corr_sm);
            FA_corr_sm_x =  [FA_corr_sm_x; FA_corr_sm];
	    MD_corr_sm = compute_MD(D_corr_sm);
            MD_corr_sm_x =  [MD_corr_sm_x; MD_corr_sm];
	    PEV_corr_sm = compute_primary_eigvec(D_corr_sm);
            PEV_corr_sm_x =  [PEV_corr_sm_x; PEV_corr_sm];
            
	    % compute new signal with scaling the bvec + DT fit
            fy_new_dwi_signal = fy_correct_signal_sm(S_corpt,g,b,L_mat);
            [D_corr_sm_fy, exitcode] = linear_vox_fit(1,fy_new_dwi_signal,g, b);
	    % FA + MD + PVE anisotropy method
            FA_corr_sm_fy = compute_FA(D_corr_sm_fy);
            FA_corr_sm_fy_x =  [FA_corr_sm_fy_x; FA_corr_sm_fy];
	    MD_corr_sm_fy = compute_MD(D_corr_sm_fy);
            MD_corr_sm_fy_x =  [MD_corr_sm_fy_x; MD_corr_sm_fy];
	    PEV_corr_sm_fy = compute_primary_eigvec(D_corr_sm_fy);
            PEV_corr_sm_fy_x =  [PEV_corr_sm_fy_x; PEV_corr_sm_fy];

	    % DT for simulate tensor without corpt
            [D_sim, exitcode] = linear_vox_fit(1,S,g, b);
            FA_sim = compute_FA(D_sim);
            FA_sim_x = [FA_sim_x; FA_sim];

            %plotDTI(DT_true);
            lala = lala + 1;
            disp(lala);

        end
    end
end
