function run_ten_sim(bvec_file,bval_file,SNR,suf_name)
    % FA/MD/PEV after and before gradient non-linear fields in Tensor simulation
    MD = 0.0008;
    n = 50;
    nv = 19406;

    % bval and bvec
    g = importdata(bvec_file);
    b = importdata(bval_file);

    % Declare array
    FA_sim_corpt_x = zeros(nv*n*n,1);
    FA_sim_noise_x = zeros(nv*n*n,1);
    FA_true_x = zeros(nv*n*n,1);
    FA_corr_bx_x = zeros(nv*n*n,1);
    FA_corr_noise_bx_x = zeros(nv*n*n,1);
    FA_corr_sm_fy_x = zeros(nv*n*n,1);
    FA_corr_noise_sm_fy_x = zeros(nv*n*n,1);

    MD_sim_corpt_x = zeros(nv*n*n,1);
    MD_sim_noise_x = zeros(nv*n*n,1);
    MD_true_x = zeros(nv*n*n,1);
    MD_corr_bx_x = zeros(nv*n*n,1);
    MD_corr_noise_bx_x = zeros(nv*n*n,1);
    MD_corr_sm_fy_x = zeros(nv*n*n,1);
    MD_corr_noise_sm_fy_x = zeros(nv*n*n,1);

    PEV_sim_corpt_x = zeros(nv*n*n,3);
    PEV_sim_noise_x = zeros(nv*n*n,3);
    PEV_true_x = zeros(nv*n*n,3);
    PEV_corr_bx_x = zeros(nv*n*n,3);
    PEV_corr_noise_bx_x = zeros(nv*n*n,3);
    PEV_corr_sm_fy_x = zeros(nv*n*n,3);
    PEV_corr_noise_sm_fy_x = zeros(nv*n*n,3);

    AD_sim_corpt_x = zeros(nv*n*n,1);
    AD_sim_noise_x = zeros(nv*n*n,1);
    AD_true_x = zeros(nv*n*n,1);
    AD_corr_bx_x = zeros(nv*n*n,1);
    AD_corr_noise_bx_x = zeros(nv*n*n,1);
    AD_corr_sm_fy_x = zeros(nv*n*n,1);
    AD_corr_noise_sm_fy_x = zeros(nv*n*n,1);

    RD_sim_corpt_x = zeros(nv*n*n,1);
    RD_sim_noise_x = zeros(nv*n*n,1);
    RD_true_x = zeros(nv*n*n,1);
    RD_corr_bx_x = zeros(nv*n*n,1);
    RD_corr_noise_bx_x = zeros(nv*n*n,1);
    RD_corr_sm_fy_x = zeros(nv*n*n,1);
    RD_corr_noise_sm_fy_x = zeros(nv*n*n,1);


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
                %Ldet_all(idx) =  det_vals(o);

                % Get L matrix at position xyz
                xyz = locations(:,o);
                L_mat(:,:) = vL(:,:,xyz(1),xyz(2),xyz(3));

                % Simulate the tensor for given PEV, FA
                x = particle_cord(1);
                y = particle_cord(2);
                z = particle_cord(3);
                DT_true = ten_compute_dt([x y z], MD,  lprep(p));

                % Simulate signal for fake tensor with L mat 
                [S_corpt, S, abvec, abval] = ten_compute_corrput_signal(DT_true, L_mat, g, b,SNR);

                % Compute FA + MD + PEV before corrpution 
                FA_true = compute_FA(DT_true);
                FA_true_x(idx) =  FA_true;
                MD_true = compute_MD(DT_true);
                MD_true_x(idx) = MD_true;
                AD_true = compute_AD(DT_true);
                AD_true_x(idx) = AD_true;
                RD_true = compute_RD(DT_true);
                RD_true_x(idx) = RD_true;
                PEV_true = compute_primary_eigvec(DT_true);
                PEV_true_x(idx,:) =  PEV_true;

                % DT fit corpt signal
                [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
                % Compute FA + MD + PEV after corrpution
                FA_sim_corpt = compute_FA(D_sim_corpt);
                FA_sim_corpt_x(idx) =   FA_sim_corpt;
                MD_sim_corpt = compute_MD(D_sim_corpt);
                MD_sim_corpt_x(idx) =   MD_sim_corpt;
                AD_sim_corpt = compute_AD(D_sim_corpt);
                AD_sim_corpt_x(idx) =   AD_sim_corpt;
                RD_sim_corpt = compute_RD(D_sim_corpt);
                RD_sim_corpt_x(idx) =   RD_sim_corpt;
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
                AD_corr_bx = compute_AD(D_corr_bx);
                AD_corr_bx_x(idx) =  AD_corr_bx;
                RD_corr_bx = compute_RD(D_corr_bx);
                RD_corr_bx_x(idx) =  RD_corr_bx;
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
                AD_corr_sm_fy = compute_AD(D_corr_sm_fy);
                AD_corr_sm_fy_x(idx) =   AD_corr_sm_fy;
                RD_corr_sm_fy = compute_RD(D_corr_sm_fy);
                RD_corr_sm_fy_x(idx) =   RD_corr_sm_fy;
                PEV_corr_sm_fy = compute_primary_eigvec(D_corr_sm_fy);
                PEV_corr_sm_fy_x(idx,:) =  PEV_corr_sm_fy;

                %noise
                % DT fit corpt signal
                [D_sim_noise, exitcode] = linear_vox_fit(1,S,g,b);
                % Compute FA + MD + PEV after corrpution
                FA_sim_noise = compute_FA(D_sim_noise);
                FA_sim_noise_x(idx) =   FA_sim_noise;
                MD_sim_noise = compute_MD(D_sim_noise);
                MD_sim_noise_x(idx) =   MD_sim_noise;
                AD_sim_noise = compute_AD(D_sim_noise);
                AD_sim_noise_x(idx) =   AD_sim_noise;
                RD_sim_noise = compute_RD(D_sim_noise);
                RD_sim_noise_x(idx) =   RD_sim_noise;
                PEV_sim_noise = compute_primary_eigvec(D_sim_noise);
                PEV_sim_noise_x(idx,:) =   PEV_sim_noise;

                % DT fit full corrected with adjected bvec and adjected bval
		[abvec_noise, abval_noise] = ten_compute_emp_noise_LR(g,b);
                [D_corr_noise_bx, exitcode] =  linear_vox_fit(1,S,abvec_noise, abval_noise);
                % compute_ADC_(D_corr_bx, g, b);
                % FA + MD + PVE full correction method
                FA_corr_noise_bx = compute_FA(D_corr_noise_bx);
                FA_corr_noise_bx_x(idx) =  FA_corr_noise_bx;
                MD_corr_noise_bx = compute_MD(D_corr_noise_bx);
                MD_corr_noise_bx_x(idx) =  MD_corr_noise_bx;
                AD_corr_noise_bx = compute_AD(D_corr_noise_bx);
                AD_corr_noise_bx_x(idx) =  AD_corr_noise_bx;
                RD_corr_noise_bx = compute_RD(D_corr_noise_bx);
                RD_corr_noise_bx_x(idx) =  RD_corr_noise_bx;
                PEV_corr_noise_bx = compute_primary_eigvec(D_corr_noise_bx);
                PEV_corr_noise_bx_x(idx,:) = PEV_corr_noise_bx;

                % Compute new signal with scaling the bvec + DT fit
                new_dwi_signal = ten_compute_noise_approx_corr(S,g,b);
                [D_corr_noise_sm_fy, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
                % FA + MD + PVE approximate correction method
                FA_corr_noise_sm_fy = compute_FA(D_corr_noise_sm_fy);
                FA_corr_noise_sm_fy_x(idx) =   FA_corr_noise_sm_fy;
                MD_corr_noise_sm_fy = compute_MD(D_corr_noise_sm_fy);
                MD_corr_noise_sm_fy_x(idx) =   MD_corr_noise_sm_fy;
                AD_corr_noise_sm_fy = compute_AD(D_corr_noise_sm_fy);
                AD_corr_noise_sm_fy_x(idx) =  AD_corr_noise_sm_fy;
                RD_corr_noise_sm_fy = compute_RD(D_corr_noise_sm_fy);
                RD_corr_noise_sm_fy_x(idx) =  RD_corr_noise_sm_fy;
                PEV_corr_noise_sm_fy = compute_primary_eigvec(D_corr_noise_sm_fy);
                PEV_corr_noise_sm_fy_x(idx,:) =  PEV_corr_noise_sm_fy;

                if mod(idx,24257500)==0
                    save("/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_"+ suf_name + "_" + int2str(idx) +".mat",'-v7.3');
                end
            end
        end
    end
    save("/nfs/masi/kanakap/projects/LR/ten_sim/CorrectionLReffectsresults_"+ suf_name  + ".mat",'-v7.3');
end


