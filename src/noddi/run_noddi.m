%addpath '/Users/kanakap/VUFall/DWI Project/code/HowToRenderVideo-master/matlab'
% addpath(genpath('../../external/spm_read_nii'))
% addpath(genpath('../..//external/spm_reslice'))
load_btable
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/nifti_matlab-master'))
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/NODDI_toolbox_v1.04'))
addpath(genpath('/nfs/masi/caily/projects/masiver/code/calc_metrics/BCT'))
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/fitting')
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/models')
addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/NODDI_toolbox_v1.05/ParforProgMonv3')

noddi_dir = '/home/local/VANDERBILT/kanakap/gradtensor/src/noddi/noddi_output_corr';
mkdir(noddi_dir);
noddi_dwi_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/sub-cIVs001_ses-s1Bx2_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.nii';
noddi_mask_file = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/prequal_dwi_cat/mrtrix_mask.nii';

noddi_data_file = fullfile(noddi_dir, 'noddi_data.mat');
CreateROI(noddi_dwi_file, noddi_mask_file, noddi_data_file);
noddi_model = MakeModel('WatsonSHStickTortIsoV_B0');
noddi_params_file = fullfile(noddi_dir, 'noddi_params.mat');

num_threads = 12;
if num_threads > 1
    batch_fitting_pk(noddi_data_file, bval_vols, bvec_vols, noddi_model, noddi_params_file, num_threads)
else
    batch_fitting_single(noddi_data_file, noddi_corrected_gradients, noddi_model, noddi_params_file)
end 

noddi_out_prefix = fullfile(noddi_dir, 'noddi');
SaveParamsAsNIfTI(noddi_params_file, noddi_data_file, noddi_dwi_file, noddi_out_prefix);
