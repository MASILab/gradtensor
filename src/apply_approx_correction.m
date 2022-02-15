function apply_approx_correction(varargin)

% Parse inputs (defaults specified here)
P = inputParser;
addOptional(P,'Limg_file','/INPUTS/L.nii.gz');
addOptional(P,'refimg_file','/INPUTS/L.nii.gz');
addOptional(P,'bval_file','/INPUTS/bval.txt');
addOptional(P,'bvec_file','/INPUTS/bvec.txt');
%addOptional(P,'mask_path','/INPUTS/mask.nii.gz');
%addOptional(P,'out_name','/INPUTS/out_name.txt');
addOptional(P,'out_dir','/OUTPUTS');
parse(P,varargin{:});

% Apply the approximate correction to signal
compute_approx_corr_dwi( ...
	P.Results.Limg_file, ...
	P.Results.refimg_file, ...
	P.Results.bval_file, ...
	P.Results.bvec_file, ...
	P.Results.mask_path, ...
	P.Resutls.out_name, ...
	P.Results.out_dir ...
	);

if isdeployed()
	exit()
end
