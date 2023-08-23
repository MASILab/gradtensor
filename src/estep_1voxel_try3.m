%% Diffusion data in bore
clear all 
close all

% subjs = { '15','30','45'}; 
subjs = {1};

bval_all = [];
for i = 1:numel(subjs)
    bval = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/scanner_bore_space2000/org.bval']);
    bval_all(:,:,i) = bval;
end
ind_b0 = find(bval <= 0.1);
ind_non_b0 = find(bval > 0.1);

bval_all_nonb0 = [];
for i = 1:numel(subjs)
    bval_all_nonb0(:,:,i) = bval_all(:,ind_non_b0,i);
end

% dwi_all = [];
% for i = 1:numel(subjs)
%     dwi = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_dwi_worldspace.nii.gz']);
%     dwi_all(:,:,:,:,i) = dwi(:,:,:,ind_non_b0);
% end
dwi_all = [];
for i = 1:numel(subjs)
    dwi = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/world_space2000_ndimagepy/dwi_worldspace.nii.gz']);
    dwi_all(:,:,:,:,i) = dwi(:,:,:,ind_non_b0);
end

% b0_all = [];
% for i = 1:numel(subjs)
%     b0 = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_dwi_worldspace.nii.gz']);
%     b0_all(:,:,:,i) = mean(b0(:,:,:,ind_b0),4);
% end
b0_all = [];
for i = 1:numel(subjs)
    b0 = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/world_space2000_ndimagepy/meanb0_worldspace.nii.gz']);
    b0_all(:,:,:,i) = b0(:,:,:,:);
end

% bvec_all = [];
% for i = 1:numel(subjs)
%     bvec = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_worldspace.bvec']);
%     bvec_all(:,:,i) = bvec(:,ind_non_b0);
% end
bvec_all = [];
for i = 1:numel(subjs)
    bvec = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/world_space2000_ndimagepy/worldspace.bvec']);
    bvec_all(:,:,i) = bvec(:,ind_non_b0);
end

% dti_all = [];
% for i = 1:numel(subjs)
%     dti = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_temp2sim_rotated/dt_temp2subj', subjs{i} ,'.nii.gz']);
%     dti_all(:,:,:,:,i) = dti(:,:,:,:);
% end
dti_all = [];
for i = 1:numel(subjs)
    dti = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/world_space2000_ndimagepy/dt.nii.gz']);
    dti_all(:,:,:,:,i) = dti(:,:,:,:);
end

fa_all = [];
for i = 1:numel(subjs)
    fa = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/world_space2000_ndimagepy/fa.nii.gz']);
    fa_all(:,:,:,i) = fa(:,:,:,:);
end

%% Simulated rotated brains
clear all 
close all


subjs = { '15','30','45'}; 

% diffusion data in bore space
bval_all = [];
for i = 1:numel(subjs)
    bval = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/scanner_bore_space2000/org.bval']);
    bval_all(:,:,i) = bval;
end
ind_b0 = find(bval <= 0.1);
ind_non_b0 = find(bval > 0.1);

bval_all_nonb0 = [];
for i = 1:numel(subjs)
    bval_all_nonb0(:,:,i) = bval_all(:,ind_non_b0,i);
end

dwi_all = [];
for i = 1:numel(subjs)
    dwi = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_dwi_worldspace.nii.gz']);
    dwi_all(:,:,:,:,i) = dwi(:,:,:,ind_non_b0);
end


b0_all = [];
for i = 1:numel(subjs)
    b0 = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_meanb0_worldspace.nii.gz']);
    b0_all(:,:,:,i) = b0(:,:,:,:);
end


bvec_all = [];
for i = 1:numel(subjs)
    bvec = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_worldspace.bvec']);
    bvec_all(:,:,i) = bvec(:,ind_non_b0);
end


dti_all = [];
for i = 1:numel(subjs)
    dti = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_temp2sim_rotated/dt_temp2subj', subjs{i} ,'.nii.gz']);
    dti_all(:,:,:,:,i) = dti(:,:,:,:);
end
% dti_all = [];
% for i = 1:numel(subjs)
%     dti = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_sim_rotated/', subjs{i} ,'_dt.nii.gz']);
%     dti_all(:,:,:,:,i) = dti(:,:,:,:);
% end

fa_all = [];
for i = 1:numel(subjs)
    fa = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/worldspace_temp2sim_rotated/fa_temp2subj', subjs{i} ,'.nii.gz']);
    fa_all(:,:,:,i) = fa(:,:,:,:);
end
i = 45; j = 62; k = 26;
%% Kids
clear all 
close all

subjs = { 'sub-cIVs001/ses-s1Bx2/','sub-cIVs002/ses-s1Bx2/','sub-cIVs006/ses-s1Bx2/', 'sub-cIVs010/ses-s1Bx3/','sub-cIVs015/ses-s1Bx2/', ...
           'sub-cIVs017/ses-s1Bx2/', 'sub-cIVs023/ses-s1Bx1/', ...
           'sub-cIVs024/ses-s1Bx1/', 'sub-cIVs025/ses-s1Bx1/','sub-cIVs029/ses-s1Bx1/','sub-cIVs031/ses-s1Bx1/', ... 
          'sub-cIVs033/ses-s1Bx2/', 'sub-cIVs035/ses-s1Bx2/','sub-cIVs037/ses-s1Bx1/', 'sub-cIVs040/ses-s1Bx1/'}%, ...
%           'sub-cIVs041/ses-s1Bx2/','sub-cIVs044/ses-s1Bx1/','sub-cIVs046/ses-s1Bx2/','sub-cIVs047/ses-s1Bx1/', ...
%           'sub-cIVs049/ses-s1Bx2/'};%,'sub-cIVs051/ses-s1Bx1/', 'sub-cIVs052/ses-s1Bx1/','sub-cIVs053/ses-s1Bx1/','sub-cIVs054/ses-s1Bx1/'};

subjs = { 'sub-cIVs006/ses-s1Bx2/','sub-cIVs024/ses-s1Bx1/', 'sub-cIVs025/ses-s1Bx1/','sub-cIVs035/ses-s1Bx2/','sub-cIVs037/ses-s1Bx1/','sub-cIVs049/ses-s1Bx2/'};      
      
% mask combined
mask_all = [];
for i = 1:numel(subjs)
    mask = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_space2000_ndimagepy/wm2dwiworld.nii.gz']);
    mask_all(:,:,:,i) = mask;
end
mask_combined = sum(mask_all, 4);

% diffusion data in bore space
bval_all = [];
for i = 1:numel(subjs)
    bval = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/scanner_bore_space2000/org.bval']);
    bval_all(:,:,i) = bval;
end
ind_b0 = find(bval <= 0.1);
ind_non_b0 = find(bval > 0.1);

bval_all_nonb0 = [];
for i = 1:numel(subjs)
    bval_all_nonb0(:,:,i) = bval_all(:,ind_non_b0,i);
end

dwi_all = [];
for i = 1:numel(subjs)
    dwi = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_space2000_ndimagepy/dwi_worldspace.nii.gz']);
    dwi_all(:,:,:,:,i) = dwi(:,:,:,ind_non_b0);
end

b0_all = [];
for i = 1:numel(subjs)
    b0 = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_space2000_ndimagepy/meanb0_worldspace.nii.gz']);
    b0_all(:,:,:,i) = b0(:,:,:,:);
end

bvec_all = [];
for i = 1:numel(subjs)
    bvec = importdata(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_space2000_ndimagepy/worldspace.bvec']);
    bvec_all(:,:,i) = bvec(:,ind_non_b0);
end

dti_all = [];
for i = 1:numel(subjs)
    dti = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_temp2subj_space2000_ndimagepy/dt_temp2subj.nii.gz']);
    dti_all(:,:,:,:,i) = dti(:,:,:,:);
end

fa_all = [];
for i = 1:numel(subjs)
    fa = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_temp2subj_space2000_ndimagepy/fa_temp2subj.nii.gz']);
    fa_all(:,:,:,i) = fa(:,:,:,:);
end

%% mask
count = 0;
for i = 1:112
    for j = 1:112
        for k = 1:54
            if mask_combined(i, j, k) == 6 
                disp([i,j,k]);
                count = count + 1;
            end
        end
    end
end
count
%% plot the signals 
close all
i = 44; j = 54; k = 32;
figure
for subj = 1:3
    data_point = [i, j, k];
    signal = dwi_all(i,j,k,:,subj);
    signal = squeeze(signal);
    plot(signal,signal,'.', 'MarkerSize', 10)
    legend()
    hold on
end
%% Plot the true tensor
% tensor back to bore space
colors = {'b','m','k','g','y'};
i = 50; j = 53; k = 27;
figure
for subj = 1:6
    disp('-----------')
    disp(subj)
    data_point = [i, j, k];
    d = dti_all(i,j,k,:,subj);
    DT_back_to_bore_true = [d(1), d(4), d(5); d(4), d(2), d(6); d(5), d(6), d(3)];
    disp(DT_back_to_bore_true)
    plotDTI_colordiff(DT_back_to_bore_true,0.2,'b')
    disp(DT_back_to_bore_true)
end

%% fa of true tensor
i = 50; j = 53; k = 27;
for subj = 1:6
    disp('-----------')
    disp(subj)
    f = fa_all(i,j,k,subj);
    disp(f)
end
%% goodness fit
for subj = 1:24
    disp('-----------')
    data_point = [i, j, k];
    d = dti_all(i,j,k,:,subj);
    DT_back_to_bore_true = [d(1), d(4), d(5); d(4), d(2), d(6); d(5), d(6), d(3)];
    
    signal = dwi_all(i,j,k,:,subj);
    signal = squeeze(signal);
    
    g = bvec_all(:,:,subj);
    b = bval_all_nonb0(:,:,subj);
    
    estimated_s = compute_signal(DT_back_to_bore_true, g, b);
    observed_counts = histcounts(estimated_s);
    expected_counts = histcounts(signal);
    [chi2stat, p_value]= chi2gof(1:length(observed_counts),'Frequency',observed_counts, 'Expected',expected_counts);
    disp(chi2stat)
    disp(p_value)
end
%% L_pred(i,j,k,:)
tic
i = 45; j = 62; k = 26;
L_pred = test_Lvol_wholeimg2(1,dti_all(i,j,k,:,:),bvec_all,bval_all_nonb0,dwi_all(i,j,k,:,:),b0_all(i,j,k,:));
disp(reshape(L_pred,3,3))
toc

%% works one voxel simulated rotated
i = 45; j = 62; k = 26;
d =  (dti_all(i,j,k,:,1));
Dxx = d(1);
Dxy = d(4);
Dxz = d(5);
Dyy = d(2);
Dyz = d(6);
Dzz = d(3);
g = bvec_all(:,:,1);
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

b0 = b0_all(i,j,k,1);
S_corpt = squeeze(dwi_all(i,j,k,:,1)) / b0; 
b = squeeze(bval_all_nonb0(:,:,1));
L0 = diag([1 1 1]);
opts.TolFun=1e-9;
opts.TolX = 1e-9;
opts.PlotFcns=@optimplotfval;
Lopt_direct = fminsearch(@(L) myLoss2(S_corpt,b0,b,qii,L),L0(:),opts);
disp(Lopt_direct)
%% works one voxel real kids 
n_subj = 6;
% i = 45; j = 62; k = 26;
i = 50; j = 53; k = 27;
qii_stack = cell(1, n_subj);
S_stack = cell(1, n_subj);
b_stack = cell(1, n_subj);
for n = 1:n_subj
    disp(n)
    d =  (dti_all(i,j,k,:,n));
    Dxx = d(1);
    Dxy = d(4);
    Dxz = d(5);
    Dyy = d(2);
    Dyz = d(6);
    Dzz = d(3);
    g = bvec_all(:,:,n);
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
    qii_stack{n} = qii;
    
    b0 = b0_all(i,j,k,n);
    S_corpt = squeeze(dwi_all(i,j,k,:,n)) / b0; 
    S_stack{n} = S_corpt;
    
    b = squeeze(bval_all_nonb0(:,:,n));
    b_stack{n} = b';
end
final_qii = vertcat(qii_stack{:});
final_Scorpt = vertcat(S_stack{:});
final_b = vertcat(b_stack{:});

L0 = diag([1 1 1]);
opts.TolFun=1e-9;
opts.TolX = 1e-9;
opts.PlotFcns=@optimplotfval;
Lopt_direct = fminsearch(@(L) myLoss2(final_Scorpt,b0,final_b',final_qii,L),L0(:),opts);
disp(reshape(Lopt_direct,3,3))
%%  wm region 
count = 0;
n_subj = 6;
for i = 1:112
    for j = 1:112
        for k = 1:54
            if mask_combined(i, j, k) == 6 
                qii_stack = cell(1, n_subj);
                S_stack = cell(1, n_subj);
                b_stack = cell(1, n_subj);
                for n = 1:n_subj
                    d =  (dti_all(i,j,k,:,n));
                    Dxx = d(1);
                    Dxy = d(4);
                    Dxz = d(5);
                    Dyy = d(2);
                    Dyz = d(6);
                    Dzz = d(3);
                    g = bvec_all(:,:,n);
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
                    qii_stack{n} = qii;

                    b0 = b0_all(i,j,k,n);
                    S_corpt = squeeze(dwi_all(i,j,k,:,n)) / b0; 
                    S_stack{n} = S_corpt;

                    b = squeeze(bval_all_nonb0(:,:,n));
                    b_stack{n} = b';
                end
                final_qii = vertcat(qii_stack{:});
                final_Scorpt = vertcat(S_stack{:});
                final_b = vertcat(b_stack{:});

                L0 = diag([1 1 1]);
                opts.TolFun=1e-9;
                opts.TolX = 1e-9;
                Lopt(i,j,k,:) = fminsearch(@(L) myLoss2(final_Scorpt,b0,final_b',final_qii,L),L0(:),opts);
                count = count + 1;
                disp(count)
                disp(reshape(Lopt,3,3))
            end
        end
    end
end     
%% Save L
Lxx = Lopt(:,:,:,1);
Lxy = Lopt(:,:,:,2);
Lxz = Lopt(:,:,:,3);
Lyx = Lopt(:,:,:,4);
Lyy = Lopt(:,:,:,5);
Lyz = Lopt(:,:,:,6);
Lzx = Lopt(:,:,:,7);
Lzy = Lopt(:,:,:,8);
Lzz = Lopt(:,:,:,9);

true_L = '/nfs/masi/kanakap/projects/estimate/L_image_space.nii';
Vref = spm_vol(true_L);
Vref = Vref(1);
Vout = rmfield(Vref,{'pinfo','private'});
Vout.dt(1) = spm_type('float32');
Vout.descrip = 'L matrix';

Limg_file = fullfile('/nfs/masi/kanakap/projects/estimate/kids_estimated_L.nii');
Vout.fname = Limg_file;
Vout.n(1) = 1;
spm_write_vol(Vout,Lxx);
Vout.n(1) = 2;
spm_write_vol(Vout,Lxy);
Vout.n(1) = 3;
spm_write_vol(Vout,Lxz);
Vout.n(1) = 4;
spm_write_vol(Vout,Lyx);
Vout.n(1) = 5;
spm_write_vol(Vout,Lyy);
Vout.n(1) = 6;
spm_write_vol(Vout,Lyz);
Vout.n(1) = 7;
spm_write_vol(Vout,Lzx);
Vout.n(1) = 8;
spm_write_vol(Vout,Lzy);
Vout.n(1) = 9;
spm_write_vol(Vout,Lzz);
%% True L 

true_L = '/nfs/masi/kanakap/projects/estimate/L_image_space.nii.gz';
% true_L = '/nfs/masi/kanakap/projects/estimate/L.nii';
true_L_vol = niftiread(true_L);
L_mat = squeeze(true_L_vol(i,j,k,:));

% Vref = spm_vol(true_L);
%%
reshape(L_pred,3,3)