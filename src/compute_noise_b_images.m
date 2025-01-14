function [adjbvec_files,adjbval_files,rLimg_file] = compute_noise_b_images( ...
	Limg_file, ...
	refimg_file, ...
	bval_file, ...
	bvec_file, ...
	out_dir ...
	)
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

% Load B values and vectors
bval = load(bval_file);
bvec = load(bvec_file);
nb = length(bval);

% Flip sign of X coordinate of bvec. I.e. FSL convention for bvecs: bvecs
% are always LAS while NII header world coords are always RAS.
%    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/FAQ#What_conventions_do_the_bvecs_use.3F
%    http://www.mrtrix.org/2016/04/15/bug-in-fsl-bvecs-handling/
disp('Flipping bvec X (LAS to RAS)')
bvec(1,:) = -bvec(1,:);
flipbvecs = true;

% Adjust bvalues and bvectors
adjbvec = nan(nv,3,nb);
adjbval = nan(nv,1,nb);

% Add complex guassian noise
initial_SNR = 100;
%mean_signal = mean(b0_vol(wm_mask))
mean_signal = 0.3333; % to do put b0
noise_std = mean_signal / initial_SNR;
randn('state',sum(100*clock));
real_noise = noise_std * randn([nv,1]);
img_noise = 1i * noise_std * randn([nv,1]);

for v = 1:nv
	
	% Most simply, the adjusted bvec is simply L * bvec. Here we are
	% operating in the image space.
	% vL = awgn([1 0 0 ; 0 1 0 ; 0 0 1], 50);
	vL = abs([1 0 0 ; 0 1 0 ; 0 0 1] + real_noise(v) + img_noise(v));
	ab = vL * bvec;
	
	% The bvecs were length 1 before adjustment, so now compute the length
	% change and adjust bvals accordingly. NOTE: adjust bval by the square
	% of the length, because the b value has a G^2 term but the vector
	% length is for G.
	len2 = sum(ab.^2);
	adjbval(v,:) = bval .* len2;
	
	% Re-normalize bvecs to length 1 to compensate for the b value
	% adjustment we just made. Skip cases where b=0.
	len = sqrt(sum(ab.^2));
	lenkeeps = len~=0;
	ab(:,lenkeeps) = ab(:,lenkeeps) ./ repmat(len(lenkeeps),3,1);
	adjbvec(v,:,:) = ab;
	
end

% Flip the bvec X back if we need to
if flipbvecs
	adjbvec(:,1,:) = -adjbvec(:,1,:);
end

% Save bvec image to file
adjbvec_files = [];
for b = 1:nb
	Vout = rmfield(VL(1),{'pinfo','private'});
	Vout.dt(1) = spm_type('float32');
	Vout.descrip = 'Adjusted bvec';
	adjbvec_files{b,1} = fullfile(out_dir,sprintf('bvec_%04d.nii',b));
	Vout.fname = adjbvec_files{b,1};
	
	for n = 1:3
		Vout.n(1) = n;
		spm_write_vol(Vout,reshape(adjbvec(:,n,b),Vout.dim));
	end
	
	system(['gzip -f ' Vout.fname]);

end

% Save bval image to file
adjbval_files = [];
for b = 1:nb
	Vout = rmfield(VL(1),{'pinfo','private'});
	Vout.dt(1) = spm_type('float32');
	Vout.descrip = 'Adjusted bval';
	adjbval_files{b,1} = fullfile(out_dir,sprintf('bval_%04d.nii',b));
	Vout.fname = adjbval_files{b,1};
	spm_write_vol(Vout,reshape(adjbval(:,b),Vout.dim));
	system(['gzip -f ' Vout.fname]);
end


