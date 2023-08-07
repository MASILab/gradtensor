% function run_ten_sim(bvec_file,bval_file,SNR,suf_name)
clear all;
bvec_file = '/home/local/VANDERBILT/kanakap/gradtensor/src/2000direction1shell.bvec';
% bvec_file = '/home/local/VANDERBILT/kanakap/gradtensor/src/30direction1shell.bvec';
% bvec_file = '/home/local/VANDERBILT/kanakap/gradtensor/src/61direction1shell.bvec';
MD = 0.0008;
n = 1;%0;
nv = 1;%0;
suf_name = 'test_metrics_det_100_dir200';
% bval and bvec # for different directions 
g = importdata(bvec_file);
g = g';
b =  repmat(1000, 1, 198); %198

% Declare array
%     FA_all = zeros(1,nv*n*n);
%     Ldet_all = zeros(nv*n*n,1);

% Generate evenly distributed FA
% FA = linspace(0.01,1,n);
FA = 0.75;
for i = 1:n
   lprep(i) = MD * (1 - FA(i)/ (3-2*FA(i)^2)^(1/2));
end

% Rank-ordered expected LR field and selected 100 by choosing the first through 100th percentiles
[L_det, vL] = vary_L_det();
n = 10;
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
%     load('theta_sphere100.mat');
%     load('phi_sphere100.mat');
load('vertices_30.mat');
idx = 0;
tic;
for i = 1%:nv
    particle_cord = vertices(i, :);
    for p = 1%:n
        for o = 1%:n 
             % saving all the variables when looped
%             idx = idx + 1;
%             FA_all(idx) =  FA(p);
%             Ldet_all(idx) =  det_vals(o);

            % Get L matrix at position xyz
            xyz = locations(:,8); %o
            L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));

            % Simulate the tensor for given PEV, FA
            x = particle_cord(1);
            y = particle_cord(2);
            z = particle_cord(3);
            DT_true = ten_compute_dt([x y z], MD,  lprep(p));

            % Simulate signal for fake tensor with L mat 
            SNR = 30;
            [S_corpt, S, abvec, abval] = ten_compute_corrput_signal(DT_true, L_mat, g, b,SNR);
            
            % save truth
            DT_true_array = [DT_true(1),DT_true(2),DT_true(3),DT_true(5),DT_true(6),DT_true(9)];
%             disp(DT_true_array)
            Ltrue(i,p,o,:) = reshape(L_mat,[],9);
            Dtrue(i,p,o,:) = DT_true_array;
            
            % initalize
            D0 = [0, 0, 0, 0, 0, 0];
            DL0 = [1, 0, 0, 0,  1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
            L0 = diag([1 1 1]);
            opts.TolFun=1e-9;
            opts.TolX = 1e-9;
%             opts.Display='iter';
%             opts.PlotFcns=@optimplotfval;
            
            % learn D first 
            Dopt = fminsearch(@(D) myLossDopt(S_corpt,b,g,D),D0(:),opts);
            Dopt_all(i,p,o,:) = Dopt;
            
            % learn L after D 
            Lopt = fminsearch(@(L) myLossDoptL(S_corpt,b,g,Dopt,L),L0(:),opts);
            Lopt_all(i,p,o,:) = Lopt;
            
            % learn them simulantous
            DLopt = fminsearch(@(DL) myLossDL(S_corpt,b,g,DL),DL0(:),opts);
            DLopt_all(i,p,o,:) = DLopt;
            
            % learn only L 
            d = DT_true;
            Dxx = d(1);
            Dxy = d(2);
            Dxz = d(3);
            Dyy = d(5);
            Dyz = d(6);
            Dzz = d(9);
            gx = g(1,:)';
            gy = g(2,:)';
            gz = g(3,:)';
            qii = [Dxx.*gx.^2  2.*Dxx.*gx.*gy  2.*Dxx.*gx.*gz  2.*Dxy.*gx.^2  2.*Dxy.*gx.*gy  2.*Dxy.*gx.*gz  2.*Dxz.*gx.^2  ...
                2.*Dxz.*gx.*gy  2.*Dxz.*gx.*gz  Dxx.*gy.^2  2.*Dxx.*gy.*gz  2.*Dxy.*gx.*gy  2.*Dxy.*gy.^2  2.*Dxy.*gy.*gz  ...
                2.*Dxz.*gx.*gy  2.*Dxz.*gy.^2  2.*Dxz.*gy.*gz  Dxx.*gz.^2  2.*Dxy.*gx.*gz  2.*Dxy.*gy.*gz  2.*Dxy.*gz.^2  ...
                2.*Dxz.*gx.*gz  2.*Dxz.*gy.*gz  2.*Dxz.*gz.^2  Dyy.*gx.^2  2.*Dyy.*gx.*gy  2.*Dyy.*gx.*gz  2.*Dyz.*gx.^2  ...
                2.*Dyz.*gx.*gy  2.*Dyz.*gx.*gz  Dyy.*gy.^2  2.*Dyy.*gy.*gz  2.*Dyz.*gx.*gy  2.*Dyz.*gy.^2  2.*Dyz.*gy.*gz  ...
                Dyy.*gz.^2  2.*Dyz.*gx.*gz  2.*Dyz.*gy.*gz  2.*Dyz.*gz.^2  Dzz.*gx.^2  2.*Dzz.*gx.*gy  2.*Dzz.*gx.*gz  ...
                Dzz.*gy.^2  2.*Dzz.*gy.*gz  Dzz.*gz.^2];
            
            Lopt_direct = fminsearch(@(L) myLoss(S_corpt,b,qii,L),L0(:),opts);
            Lopt_direct_all(i,p,o,:) = Lopt_direct;
            
            % compute error
%             Derror(i,p,o) = norm(Dopt' - DT_true_array);
%             Lerror(i,p,o) = norm(Lopt' - reshape(L_mat,[],9));

            DerrorTogether(i,p,o) = norm(DLopt(10:end)' - DT_true_array);
            LerrorTogether(i,p,o) = norm(DLopt(1:9)' - reshape(L_mat,[],9));
            
            LerrorDirect(i,p,o) = norm(Lopt_direct' - reshape(L_mat,[],9));
            
            % save the determinant 
            det_DL(i,p,o) =  det(L_mat) - det(reshape(DLopt(1:9),3,3)) ;
            det_LDirect(i,p,o) =  det(L_mat) - det(reshape(Lopt_direct,3,3));
            
            % corpt 
            [D_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
%             disp(D_corpt)
            FA_corpt = compute_FA(D_corpt);
            FA_corpt_x(i,p,o) =   FA_corpt;
            MD_corpt = compute_MD(D_corpt);
            MD_corpt_x(i,p,o) =   MD_corpt;
            PEV_corpt = compute_primary_eigvec(D_corpt);
            PEV_corpt_x(i,p,o,:) =   PEV_corpt;
            
            % correct using predicted L from learning both D and L
            L_pred = reshape(DLopt(1:9),3,3);
            [abvec_pred, abval_pred] = ten_compute_b( g, b, L_pred);
            [D_pred_corr, exitcode] = linear_vox_fit(1,S_corpt,abvec_pred, abval_pred);
            Dunknow_pred_corr = D_pred_corr;
%             disp(D_pred_corr)
            FA_pred_corr = compute_FA(D_pred_corr);
            FA_pred_corr_x(i,p,o) =   FA_pred_corr;
            MD_pred_corr = compute_MD(D_pred_corr);
            MD_pred_corr_x(i,p,o) =   MD_pred_corr;
            PEV_pred_corr = compute_primary_eigvec(D_pred_corr);
            PEV_pred_corr_x(i,p,o,:) =   PEV_pred_corr;
            
            % correct using predicted L from learning only L
            L_pred = reshape(Lopt_direct,3,3);
            [abvec_pred, abval_pred] = ten_compute_b( g, b, L_pred);
            [D_pred_corr, exitcode] = linear_vox_fit(1,S_corpt,abvec_pred, abval_pred);
            D_Ldpred_corr = D_pred_corr;
%             disp(D_Ldpred_corr)
            FA_pred_corr = compute_FA(D_pred_corr);
            FA_Ldpred_corr_x(i,p,o) =   FA_pred_corr;
            MD_pred_corr = compute_MD(D_pred_corr);
            MD_Ldpred_corr_x(i,p,o) =   MD_pred_corr;
            PEV_pred_corr = compute_primary_eigvec(D_pred_corr);
            PEV_Ldpred_corr_x(i,p,o,:) =   PEV_pred_corr;
            
            % correct using true D
            [abvec_true, abval_true] = ten_compute_b( g, b, L_mat);
            [D_true_corr, exitcode] = linear_vox_fit(1,S_corpt,abvec_true, abval_true);
            FA_true_corr = compute_FA(D_true_corr);
            FA_true_corr_x(i,p,o) =   FA_true_corr;
            MD_true_corr = compute_MD(D_true_corr);
            MD_true_corr_x(i,p,o) =   MD_true_corr;
            PEV_true_corr = compute_primary_eigvec(D_true_corr);
            PEV_true_corr_x(i,p,o,:) =   PEV_true_corr;
           
        end
    end
end
toc
% save("/nfs/masi/kanakap/projects/estimate/simulation/results_"+ suf_name  + ".mat",'-v7.3');
 