%addpath '/Users/kanakap/VUFall/DWI Project/code/HowToRenderVideo-master/matlab'
% addpath(genpath('../../external/spm_read_nii'))
% addpath(genpath('../..//external/spm_reslice'))
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/nifti_matlab-master'))
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/NODDI_toolbox_v1.04'))
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/BCT'))
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/fitting')
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/models')
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/ParforProgMonv3')

% after merge_shells
noddi_mask_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/mrtrix_mask.nii';
gunzip('/home-nfs2/local/VANDERBILT/kanakap/gradtensor/src/noddi/noddi_output/b1k2k.nii.gz')
noddi_dwi_file = '/home-nfs2/local/VANDERBILT/kanakap/gradtensor/src/noddi/noddi_output/b1k2k.nii';
noddi_dir = '/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/noddi_output';
mkdir(noddi_dir);


noddi_gradients = FSL2Protocol(b1k2k_bvals_file, b1k2k_bvecs_file);
batch_fitting(noddi_data_file, noddi_gradients, noddi_model, noddi_params_file, num_threads)

noddi_out_prefix = fullfile(noddi_dir, 'noddi');
SaveParamsAsNIfTI(noddi_params_file, noddi_data_file, noddi_dwi_file, noddi_out_prefix);