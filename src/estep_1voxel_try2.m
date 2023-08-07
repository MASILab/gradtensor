clear all 
close all
% 'sub-cIVs005/ses-s1Bx3/'
% subjs = {'sub-cIVs001/ses-s1Bx2', 'sub-cIVs002/ses-s1Bx2', 'sub-cIVs006/ses-s1Bx2', 'sub-cIVs015/ses-s1Bx2', 'sub-cIVs017/ses-s1Bx2'};
subjs = { 'sub-cIVs001/ses-s1Bx2/','sub-cIVs002/ses-s1Bx2/','sub-cIVs006/ses-s1Bx2/', 'sub-cIVs010/ses-s1Bx3/','sub-cIVs015/ses-s1Bx2/', ...
           'sub-cIVs017/ses-s1Bx2/', 'sub-cIVs023/ses-s1Bx1/', ...
           'sub-cIVs024/ses-s1Bx1/', 'sub-cIVs025/ses-s1Bx1/','sub-cIVs029/ses-s1Bx1/','sub-cIVs031/ses-s1Bx1/', ... 
          'sub-cIVs033/ses-s1Bx2/', 'sub-cIVs035/ses-s1Bx2/','sub-cIVs037/ses-s1Bx1/', 'sub-cIVs040/ses-s1Bx1/', ...
          'sub-cIVs041/ses-s1Bx2/','sub-cIVs044/ses-s1Bx1/','sub-cIVs046/ses-s1Bx2/','sub-cIVs047/ses-s1Bx1/', ...
          'sub-cIVs049/ses-s1Bx2/'};%,'sub-cIVs051/ses-s1Bx1/', 'sub-cIVs052/ses-s1Bx1/','sub-cIVs053/ses-s1Bx1/','sub-cIVs054/ses-s1Bx1/'};
% mask combined
mask_all = [];
for i = 1:numel(subjs)
    mask = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/worldspace_fast_seg/wm2dwi.nii.gz']);
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
    b0 = niftiread(['/nfs/masi/kanakap/projects/estimate/MASiVar_kids/', subjs{i}, '/world_space2000_ndimagepy/dwi_worldspace.nii.gz']);
    b0_all(:,:,:,i) = mean(b0(:,:,:,ind_b0),4);
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

template_dti = niftiread('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs005/ses-s1Bx2/scanner_bore_space2000/dt_resamp_masked.nii.gz');

%% mask
count = 0
for i = 1:112
    for j = 1:112
        for k = 1:54
            if mask_combined(i, j, k) == 2 
                disp([i,j,k]);
                count = count + 1;
            end
        end
    end
end
count
%% plot the signals 
close all
i = 72;j = 41; k = 27;
figure
for subj = 1:10
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
% i = 56;j = 64; k = 29;
i = 44; j = 62; k = 28;
figure
for subj = 1:21
    disp('-----------')
    disp(subj)
    data_point = [i, j, k];
    d = dti_all(i,j,k,:,subj);
    DT_back_to_bore_true = [d(1), d(4), d(5); d(4), d(2), d(6); d(5), d(6), d(3)];
    disp(DT_back_to_bore_true)
    plotDTI_colordiff(DT_back_to_bore_true,0.2,'b')
    disp(DT_back_to_bore_true)
end

%%
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
% i = 56;j = 64; k = 29;
i = 44; j = 62; k = 28;
L_pred = test_Lvol_wholeimg2(20,dti_all(i,j,k,:,:),bvec_all,bval_all_nonb0,dwi_all(i,j,k,:,:),b0_all(i,j,k,:));
toc
%% True L 

true_L = '/nfs/masi/kanakap/projects/estimate/L_image_space.nii.gz';
% true_L = '/nfs/masi/kanakap/projects/estimate/L.nii';
true_L_vol = niftiread(true_L);
L_mat = squeeze(true_L_vol(i,j,k,:));

% Vref = spm_vol(true_L);

