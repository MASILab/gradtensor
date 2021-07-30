%function sim_dt_FA_MD(MD,FA)
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))


flags = struct( ...
        'mask',true, ...
        'mean',false, ...
        'interp',1, ...
        'which',1, ...
        'wrap',[0 0 0], ...
        'prefix','r' ...
        );

m = 100;
FA_sim_corpt_x = zeros(size(1:m));
FA_true_x = zeros(size(1:m));
FA_sim_x = zeros(size(1:m));
FA_corr_bx_x = zeros(size(1:m));
FA_corr_sm_x = zeros(size(1:m));

%Limg_file = '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz';
%refimg_file ='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot_Lest_sig.nii';
out_dir ='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap';

mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii';
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);
% Unzip
% if strcmp(Limg_file(end-2:end),'.gz')
%         system(['gunzip -kf ' Limg_file]);
%         Limg_file = Limg_file(1:end-3);
% end
% if strcmp(refimg_file(end-2:end),'.gz')
%         system(['gunzip -kf ' refimg_file]);
%         refimg_file = refimg_file(1:end-3);
% end
%spm_reslice({refimg_file; Limg_file},flags);
%[p,n,e] = fileparts(Limg_file);
%rLimg_file = fullfile(p,['r' n e]);
%movefile(rLimg_file,fullfile(out_dir,'L_resamp.nii'));
rLimg_file = fullfile(out_dir,'L_resamp.nii');
VL = spm_vol(rLimg_file);
L = spm_read_vols(VL);
%L = reshape(L,[],9);
nv = size(L,1);
vL = zeros(3,3,size(L,1),size(L,2),size(L,3));
vL(1,1,:,:,:) = L(:,:,:,1);
vL(1,2,:,:,:) = L(:,:,:,2);
vL(1,3,:,:,:) = L(:,:,:,3);
vL(2,1,:,:,:) = L(:,:,:,4);
vL(2,2,:,:,:) = L(:,:,:,5);
vL(2,3,:,:,:) = L(:,:,:,6);
vL(3,1,:,:,:) = L(:,:,:,7);
vL(3,2,:,:,:) = L(:,:,:,8);
vL(3,3,:,:,:) = L(:,:,:,9);
%L_det = zeros(size(L,1),size(L,2),size(L,3));
%L_det = reshape(L_det,[],96*96*68);
%L_det = zeros(size(1:96*96*68));
for x = 1:size(L,1)
        for y = 1:size(L,2)
            for z = 1:size(L,3)
                if mask_vol(x,y,z)
                    %dwmri_vols(x,y,z,i) = rot90(dwmri_vols(x,y,z,i));
                    L_mat = squeeze(vL(:,:,x,y,z));
                    if det(L_mat(:,:)) ~= 0
                        L_det(x,y,z) = det(L_mat(:,:));
                    end
                end
            end
        end
end
%end
%{
maxi = 96;
i = randi([1,maxi],1,m);

maxj = 96;
j = randi([1,maxj],1,m);

maxk = 68;
k = randi([1,maxk],1,m);

for v = 1:m
    
    L_mat = vL(:,:,i(v),j(v),k(v));
    MD = 0.0008;
    

    %MD = 0.0008;
    % Trace fixed
    %TrD = 0.0021; %derek jones
    %MD = TrD/3;
    
    % FA change
    FA = 0.3;
    %l1 = (MD) * (1 + (2*FA) / sqrt( 3 - (2 * FA^2)));
    %lprep = (MD) * (1 - FA / sqrt( 3 - (2 * FA^2)));
    lprep = MD * (1 - FA/ (3-2*FA^2)^(1/2));
% for PVE
    %for i = 1:n
     
    %a1 = 0.0;
    %b1 = 0.1;
    %sz = [1 n];

    theta = 60; %random 0-180 2 * pi * unifrnd(a1,b1,sz);
    phi = 30; % 0-90 acos(1 - 2 * unifrnd(a1,b1,sz));
        
 
        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);
        PVE = [x y z];
        
        DT_true = sim_dt(PVE, MD,  lprep);

        
         [S_corpt, S, abvec, abval, g, b] = sim_signal(DT_true, L_mat);
         
         FA_true = compute_FA(DT_true);
         FA_true_x(v) = FA_true;
         
        [D_sim_corpt, exitcode] = linear_vox_fit(1,S_corpt,g, b);
        FA_sim_corpt = compute_FA(D_sim_corpt);
        FA_sim_corpt_x(v) =  FA_sim_corpt;
        
        [D_corr_bx, exitcode] =  linear_vox_fit(1,S_corpt,abvec, abval);
        FA_corr_bx = compute_FA(D_corr_bx);
        disp(FA_corr_bx)
        FA_corr_bx_x(v) = FA_corr_bx;

        new_dwi_signal = correct_signal_sm(S_corpt,b);
        [D_corr_sm, exitcode] = linear_vox_fit(1,new_dwi_signal,g, b);
        FA_corr_sm = compute_FA(D_corr_sm);
        FA_corr_sm_x(v) =  FA_corr_sm;

        [D_sim, exitcode] = linear_vox_fit(1,S,g, b);
        FA_sim = compute_FA(D_sim);
        FA_sim_x(v) = FA_sim;        
        
        %FA = compute_FA(D);
        %plotDTI(DT_true,0.002);
end
%end
%}
