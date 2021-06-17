% convert bvec and bvec to volumes
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')

% Load data
dwi_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.nii';
dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
b0_vol = dwi_vols(:,:,:,1);
dwi_vols = dwi_vols(:,:,:,1:end);
    
% get the size from dwmri
VD = spm_vol(dwi_path);
D = spm_read_vols(VD);
disp(size(D));
D = reshape(D,[],25);
nv = size(D,1);

% load 
bvec_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.bvec';
bval_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.bval';
bvec = load(bvec_path);
bvec = bvec()
bval = load(bval_path);
nb = length(bval);
disp(size(bvec))
disp(size(bval))

%{
% for bvec
out_dir = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/bvecs';
adjbvec = nan(nv,3,nb);
bvec_files = [];
for b = 1:nb
    Vout = rmfield(VD(1),{'pinfo','private'});
	Vout.dt(1) = spm_type('float32');
	Vout.descrip = 'dwmri bvec vol';
    bvec_files{b,1} = fullfile(out_dir,sprintf('bvec_%04d.nii',b));
    Vout.fname = bvec_files{b,1};
    
    for n = 1:3
		Vout.n(1) = n;
		spm_write_vol(Vout,reshape(adjbvec(:,n,b),Vout.dim));
    end
	
	system(['gzip -f ' Vout.fname]);
    
end


out_dir = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/bvals';
adjbval = nan(nv,1,nb);
bval_files = [];
for b = 1:nb
    Vout = rmfield(VD(1),{'pinfo','private'});
	Vout.dt(1) = spm_type('float32');
	Vout.descrip = 'dwmri bval vol';
    bval_files{b,1} = fullfile(out_dir,sprintf('bval_%04d.nii',b));
    Vout.fname = bval_files{b,1};
    spm_write_vol(Vout,reshape(adjbval(:,b),Vout.dim));
	system(['gzip -f ' Vout.fname]);
    
end
%}

%bvecs = [];
%bvecs = cat(4,bvecs,permute(nifti_utils.load_untouch_nii_vol_scaled(bvec_path,'double'),[1 2 3 5 4]));



%bvals = [];
%bvals = cat(4,bvals,nifti_utils.load_untouch_nii_vol_scaled(bval_path,'double'));


mask_vol = '';

% Initialize DT and exitcode volumes
exitcode_vol = zeros(size(b0_vol));
eig_vol = zeros(size(b0_vol,1),size(b0_vol,2),size(b0_vol,3),3);
primary_vec_vol = zeros(size(eig_vol)); 

[DT_mat, exitcode] = linear_vol_fit(b0_vol,dwi_vols,bvec,bval,mask_vol);

if sum(isnan(DT_mat(:))) == 0 && sum(isinf(DT_mat(:))) == 0
    [v, e] = eig(DT_mat);
    e = diag(e);
    [max_eig, pos] = max(e);
    primary = v(:,pos);
						
    eig_vol(:,:,:,:) = e;
	primary_vec_vol(:,:,:,:) = primary;
    exitcode_vol(:,:,:) = exitcode;
end

out_dir = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/';
out_name = 'try';
MD = (eig_vol(:,:,:,1) + eig_vol(:,:,:,2) + eig_vol(:,:,:,3))./3;
FA = sqrt(1/2) .* (sqrt( (eig_vol(:,:,:,1) - eig_vol(:,:,:,2)).^2 + (eig_vol(:,:,:,2) - eig_vol(:,:,:,3)).^2  + (eig_vol(:,:,:,3) - eig_vol(:,:,:,1)).^2 ) ./ sqrt(eig_vol(:,:,:,1).^2 + eig_vol(:,:,:,2).^2 + eig_vol(:,:,:,3).^2));
nii = load_untouch_nii(dwi_path);
nii.img = MD;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_md']),nii,'double');
nii.img = FA;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_fa']),nii,'double');
nii.img = primary_vec_vol;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_primary_eigvec']),nii,'double');
    

