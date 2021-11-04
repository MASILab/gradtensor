function bl_compute_simple_corr_dwi( ...
	Limg_file, ...
	refimg_file, ...
	bval_file, ...
	bvec_file, ...
	out_dir ...
	)

% compute the Frank's simple corrected diffusion method
% COMPUTE_B_IMAGES  Adjust B values for gradient coil nonlinearity
%
% Inputs
%    L_file           L matrix file
%    bval_file        B value text file
%    bvec_file        B vector text file
%
% Outputs
%    <refimg>_bval_*.nii    B value image for each diffusion direction
%    <refimg>_bvec_*.nii    B vector image for each diffusion direction

% Unzip
if strcmp(Limg_file(end-2:end),'.gz')
	system(['gunzip -kf ' Limg_file]); 
	Limg_file = Limg_file(1:end-3);
end
if strcmp(refimg_file(end-2:end),'.gz')
	system(['gunzip -kf ' refimg_file]); 
	refimg_file = refimg_file(1:end-3);
end

% Resample grad tensor to DTI image space
flags = struct( ...
	'mask',true, ...
	'mean',false, ...
	'interp',1, ...
	'which',1, ...
	'wrap',[0 0 0], ...
	'prefix','r' ...
	);
spm_reslice({refimg_file; Limg_file},flags);
[p,n,e] = fileparts(Limg_file);
rLimg_file = fullfile(p,['r' n e]);
movefile(rLimg_file,fullfile(out_dir,'L_resamp.nii'));
rLimg_file = fullfile(out_dir,'L_resamp.nii');

% Load the grad tensor and reshape. Initial dimensions are x,y,z,e where e
% is Lxx, Lxy, Lxz, Lyx, Lyy, etc. We need to reshape to i,j,v where Lij is
% the tensor for voxel v.
VL = spm_vol(rLimg_file);
L = spm_read_vols(VL);
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
vL = reshape(vL,[3,3,96,96,68]);


% Load data
dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(refimg_file,'double');
size(dwmri_vols)
b0_vol = dwmri_vols(:,:,:,1);
dwi_vols = dwmri_vols(:,:,:,2:end);
size(dwi_vols)

mask_path = '/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/mask.nii';
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
mask_vol = logical(mask_vol);

% Load B values and vectors
bval = load(bval_file);
bvec = load(bvec_file);
bval = bval(2:end);
bvec = bvec(1:3,2:end); %remove b0
nb = length(bval);


system(['gzip -f ' rLimg_file]);

new_dwi_signal = zeros(size(dwi_vols));
size(new_dwi_signal)

for i = 1:24
    for x = 1:size(dwi_vols,1)
        for y = 1:size(dwi_vols,2)
            for z = 1:size(dwi_vols,3)
                if mask_vol(x,y,z)
			% LR 
                        L_mat = squeeze(vL(:,:,x,y,z));
			og = bvec(:,i);
                        ob = bval(i);
                        og(1) = -og(1);

			% adjus bvec by L*bvec and then compute the length change 
                        gg = L_mat * og;
                        norm_gg = sum(gg.^2);
                        
			% compute the new signal with length change 
			new_dwi_signal(x,y,z,i) = (b0_vol(x,y,z)) * (exp( (log(dwi_vols(x,y,z,i)/b0_vol(x,y,z))) / norm_gg));
                        ADC = log(new_dwi_signal(x,y,z,i) / (b0_vol(x,y,z))) * (1 / (b0_vol(x,y,z) - ob));
                        fprintf('ADC %f for volume %i\n', [ADC, i]);
                end
            end
        end
    end
end

size(new_dwi_signal)
out_name = 'Lest_imagespace';
corrected_signal = zeros(size(dwmri_vols));
corrected_signal(:,:,:,1) = b0_vol ;
corrected_signal(:,:,:,2:end) = new_dwi_signal ;
nii = load_untouch_nii(refimg_file);
nii.img = corrected_signal;
%nii.img = corrected_signal;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_corrected_sig']),nii,'double');

