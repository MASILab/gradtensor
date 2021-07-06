% fit dwmri to dti
% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')


dwi_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.nii';
dwi_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_path,'double');
b0_vol = dwi_vols(:,:,:,1);
dwi_vols = dwi_vols(:,:,:,2:end);

VD = spm_vol(dwi_path);
D = spm_read_vols(VD);
disp(size(D));
D = reshape(D,[],25);
nv = size(D,1);

bvec_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bvec';
bval_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bval';
bvec = load(bvec_path);
bval = load(bval_path);
nb = length(bval);

out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/bvecs';
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


out_dir = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/bvals';
adjbval = nan(nv,1,nb);
bvec_files = [];
for b = 1:nb
    Vout = rmfield(VD(1),{'pinfo','private'});
	Vout.dt(1) = spm_type('float32');
	Vout.descrip = 'dwmri bval vol';
    bval_files{b,1} = fullfile(out_dir,sprintf('bval_%04d.nii',b));
    Vout.fname = bval_files{b,1};
    
    for n = 1:3
		Vout.n(1) = n;
		spm_write_vol(Vout,reshape(adjbval(:,n,b),Vout.dim));
    end
	
	system(['gzip -f ' Vout.fname]);
    
end


%bvecs = [];
%bvecs = cat(4,bvecs,permute(nifti_utils.load_untouch_nii_vol_scaled(bvec_path,'double'),[1 2 3 5 4]));
%disp(size(bvecs))


%bvals = [];
%bvals = cat(4,bvals,nifti_utils.load_untouch_nii_vol_scaled(bval_path,'double'));
%disp(size(bvals))

mask_vol = '';

%linear_vol_fit(b0_vol,dwi_vols,bvecs,bvals,mask_vol);
