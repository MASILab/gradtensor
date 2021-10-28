%function sim_dt_FA_MD(MD,FA)
function [L_det, vL] = vary_L_det()
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')
% for spm
addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))



%Limg_file = '/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz';
%refimg_file ='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot_Lest_sig.nii';
out_dir ='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap';

mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/mask.nii';
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);

rLimg_file = fullfile(out_dir,'L_resamp.nii');
VL = spm_vol(rLimg_file);
L = spm_read_vols(VL);
size_x = size(L,1);
size_y = size(L,2);
size_z = size(L,3);
L = reshape(L,[],9);
nv = size(L,1);
vL = zeros(3,3,nv);
vL(1,1,:) = L(:,1);
vL(1,2,:) = L(:,2);
vL(1,3,:) = L(:,3);
vL(2,1,:) = L(:,4);
vL(2,2,:) = L(:,5);
vL(2,3,:) = L(:,6);
vL(3,1,:) = L(:,7);
vL(3,2,:) = L(:,8);
vL(3,3,:) = L(:,9);
vL = reshape(vL,[3,3,size_x,size_y,size_z]);
L_det = zeros(size_x,size_y,size_z);
%L_det = reshape(L_det,[],96*96*68);
%L_det = zeros(size(1:96*96*68));
for x = 1:size_x
        for y = 1:size_y
            for z = 1:size_z
                %if mask_vol(x,y,z)
                    %dwmri_vols(x,y,z,i) = rot90(dwmri_vols(x,y,z,i));
                    L_mat = squeeze(vL(:,:,x,y,z));
                    L_det(x,y,z) = det(L_mat(:,:));
                %end
            end
        end
end
%rL_det = reshape(L_det,[],size(L_det,1)*size(L_det,2)*size(L_det,3));
%rL_det(rL_det == 0) = [];


%A= L_det;
%figure;imagesc(A(:,:,50))
%figure;imagesc(A(:,:,32))
%figure;imagesc(A(:,:,12))
%blah = A(A>0);
%[out,idx] = sort(blah,'ascend');
%indices_to_chose = round(linspace(1,119862,100));
%figure; plot(out(indices_to_chose));
%vec = out(indices_to_chose);
%for i = 1:100
%  [r,c,v] = ind2sub(size(A),find(A == vec(i)));
%  det_vals(i) = A(r,c,v);
%  locations(:,i) = [r c v];
%end
    
   

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
end
%end
%}
