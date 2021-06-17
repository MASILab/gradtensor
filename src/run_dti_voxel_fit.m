% run dti voxel fit on the corrected bvec and bvals phantoms

% for nifti_utils
addpath('/nfs/masi/kanakap/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')
% for load_untouch_header_only
addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')

gunzip('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii.gz')

dwi_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/INPUTS/dwmri.nii';
bvec_folder = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/';
bval_folder = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/';
mask_path = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii';
out_dir = '/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/';
out_name = 'p_3tb_posA_mask';

dti_voxel_fit(dwi_path,bvec_folder,bval_folder,mask_path, out_dir, out_name);
