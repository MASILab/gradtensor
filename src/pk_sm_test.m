addpath(genpath('../external/spm_read_nii'))
addpath(genpath('../external/spm_reslice'))
addpath('sh_basis')

if 0
	fieldmaps_to_gradtensor( ...
		'fieldmap_hz_0_file','/nfs/masi/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_0_Hz.nii', ...
		'fieldmap_hz_x_file','/nfs/masi/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_X_Hz.nii', ...
		'fieldmap_hz_y_file','/nfs/masi/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_Y_Hz.nii', ...
		'fieldmap_hz_z_file','/nfs/masi/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/INPUTS/B0_Z_Hz.nii', ...
		'sh_order',7, ...
		'symmetry','sym', ...
		'coord_geom','Philips_headfirst_supine', ...
		'image_radius',135, ...
		'out_dir','/nfs/masi/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS' ...
		);
end


%create output dir
apply_simple_correction( ...	
    'Limg_file','/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz',   ...
    'refimg_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ten_est/rot_Lest_sig.nii', ...	
    'bval_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bval', ...
	'bvec_file','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/INPUTS/dwmri.bvec', ...
	'out_dir','/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_simple_method_future_fieldmap' ...
	);