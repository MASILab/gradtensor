function compute_noise_approximate_corr_dwi( ...
	refimg_file, ...
	bval_file, ...
	bvec_file, ...
	mask_path, ...
	out_name, ...
	out_dir ...
	)

% compute_approximate_corr_dwi Compute the approximate corrected signal for gradient coil nonlinearity
%
% Inputs
%    refimg_file
%    bval_file        B value text file
%    bvec_file        B vector text file
%    mask_path        Path to the mask
%    out_name         Prefix of the output filename
%    out_dir          Output directory


% Load data
dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(refimg_file,'double');
bvec = load(bvec_file);
bval = load(bval_file);

% dwi and non-dwi
ind_b0 = find(~bval);
ind_non_b0 = find(bval);

all_b0_vol = dwmri_vols(:,:,:,ind_b0);
b0_vol = mean(all_b0_vol,4);
dwi_vols = dwmri_vols(:,:,:,ind_non_b0);

% Load the mask
if ~exist('mask_path','var') || isempty(mask_path)
        disp('making mask')
        mask_vol = true(size(b0_vol));
    else
        %gunzip(mask_path)
        mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_path,'double');
        mask_vol = logical(mask_vol);
    end

% Load B values and vectors
bval = bval(ind_non_b0);
bvec = bvec(:,ind_non_b0); %remove b0
nb = length(bval);

% Resmaple L image

% Create an zeros array fo new signal
new_dwi_signal = zeros(size(dwi_vols));
size(new_dwi_signal)


% Add complex guassian noise
initial_SNR = 100;
%mean_signal = mean(b0_vol(wm_mask))
mean_signal = 0.3333;
noise_std = mean_signal / initial_SNR;
randn('state',sum(100*clock));
real_noise = noise_std * randn(size(dwi_vols));
img_noise = 1i * noise_std * randn(size(dwi_vols));

% Compute the new signal by looping thought the number of bvals and along the x, y, z axis
for i = 1:size(dwi_vols,4)
    for x = 1:size(dwi_vols,1)
        for y = 1:size(dwi_vols,2)
            for z = 1:size(dwi_vols,3)
                if mask_vol(x,y,z)
                        % LR 
                        %L_mat = awgn([1 0 0 ; 0 1 0 ; 0 0 1], 50);
			L_mat = abs([1 0 0 ; 0 1 0 ; 0 0 1] + real_noise(x,y,z,i) + img_noise(x,y,z,i));
                        % Original bval and bvec
                        og = bvec(:,i);
                        ob = bval(i);
                        og(1) = -og(1);

                        % Adjus bvec by L*bvec and then compute the length change 
                        gg = L_mat * og;
                        norm_gg = sum(gg.^2);
            
                        % Compute the new signal with length change 
                        %new_dwi_signal(v) = b0 * exp( (log(S_corpt(v)/b0)) / len2 );
                        new_dwi_signal(x,y,z,i) = b0_vol(x,y,z) * exp( (log(dwi_vols(x,y,z,i)/b0_vol(x,y,z))) / norm_gg);
                        % Compute the ADC if required
                        %bv_b0 = 0;
                        %ADC = log(new_dwi_signal(x,y,z,i) / (b0_vol(x,y,z))) * (1 / (bv_b0 - ob));
                        %fprintf('ADC simple corr %f for volume %i\n', [ADC, i]);
                end
            end
        end
    end
end


% Output prefix file name
%out_name = 'Lest';
% Save the corrected signal (b0 was not used during computation so add that seperately to volume 1)
corrected_signal = zeros(size(dwmri_vols));
corrected_signal(:,:,:,ind_b0) = all_b0_vol ;
corrected_signal(:,:,:,ind_non_b0) = new_dwi_signal ;
nii = load_untouch_nii(refimg_file);
nii.img = corrected_signal;
nifti_utils.save_untouch_nii_using_scaled_img_info(fullfile(out_dir, [out_name '_sig']),nii,'double');
