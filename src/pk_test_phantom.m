addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))
addpath('sh_basis')

if 1
	fieldmaps_to_gradtensor( ...
		'fieldmap_hz_0_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_0_Hz.nii', ...
		'fieldmap_hz_x_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_X_Hz.nii', ...
		'fieldmap_hz_y_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_Y_Hz.nii', ...
		'fieldmap_hz_z_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_Z_Hz.nii', ...
		'sh_order',7, ...
		'symmetry','sym', ...
		'coord_geom','Philips_headfirst_supine', ...
		'image_radius',135, ...
		'out_dir','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/Jul27_OUTPUTS' ...
		);
end


%create output dir
%apply_gradtensor_to_b( ...	
 %   'Limg_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz',   ...
 %   'refimg_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_signal_estimate/try1_inv_Lest_sig.nii', ...	
 %   'bval_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bval', ...
%	'bvec_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.bvec', ...
%	'out_dir','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUT_est_corrected' ...
%	);
