# lr pipeline

# Inputs
#dwi_path=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.nii
dwi_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.nii
bvec_folder=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/
bval_folder=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/
org_bvec_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bvec
org_bval_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bval
mask_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii

out_dir=/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_32_1
#rL_path=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii

mkdir $out_dir
uncorr_out_name=uncorrected
initial_SNR=100
uncorrected_sig=$out_dir/uncorrected_Lest_sig.nii

Limg_file=/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz
refimg_file=$uncorrected_sig
approx_out_name=approx_corrected

approx_corrected=$out_dir/approx_corrected_sig.nii

# Reconstruct corrput signal
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_read_nii'));addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_reslice'));addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/sh_basis');addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');disp('Reconstruct corrput signal');compute_noise_corrput_signal('$dwi_path','$bvec_folder','$bval_folder','$mask_path','$out_dir','$uncorr_out_name','$Limg_file','$org_bvec_path','$org_bval_path',$initial_SNR);disp('Compute FA/MD/PEV of corrput signal');dti_voxel_fit_sig('$uncorrected_sig','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$uncorr_out_name');disp('Approximate correction');compute_approx_corr_dwi('$Limg_file','$refimg_file','$org_bval_path','$org_bvec_path','$mask_path','$approx_out_name','$out_dir'); disp('Compute FA/MD/PEV of approximate correction');dti_voxel_fit_sig('$approx_corrected','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$approx_out_name'); disp('Empirical correction');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$refimg_file','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$out_dir');exit"

emp_bvec_folder=$out_dir/emp_corrected_bvec
emp_bval_folder=$out_dir/emp_corrected_bval
mkdir $emp_bvec_folder
mkdir $emp_bval_folder
mv $out_dir/bvec_00* $emp_bvec_folder
mv $out_dir/bval_00* $emp_bval_folder
emp_out_name=emp_corrected
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('../external/spm_read_nii'));addpath(genpath('../external/spm_reslice'));addpath('sh_basis');addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');disp('Compute FA/MD/PEV of empirical correction');dti_voxel_fit_masiver('$uncorrected_sig','$emp_bvec_folder','$emp_bval_folder','$mask_path','$out_dir','$emp_out_name','$org_bval_path','$org_bvec_path');exit"

<<com
# Compute FA/MD/PEV of corrput signal
uncorrected_sig='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUT1/uncorrected_Lest_sig.nii'
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "dti_voxel_fit_sig($uncorrected_sig,$bvec_path,$bval_path,$mask_path,$out_dir,$out_name)"

# Approximate correction
Limg_file='/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz'
refimg_file=uncorrected_sig
out_name=approx_corrected
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "apply_approx_correction($Limg_file,$refimg_file,$org_bval_path,$org_bvec_path,$mask_path,$out_name,$out_dir)"

# Compute FA/MD/PEV of approximate correction
approx_corrected='/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUT1/approx_corrected_sig.nii'
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "dti_voxel_fit_sig($approx_corrected,$org_bvec_path,$org_bval_path,$mask_path,$out_dir,$out_name)"

# Empirical correction 
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "apply_gradtensor_to_b($Limg_file,$refimg_file,$org_bval_path,$org_bvec_path,$out_dir)"
emp_bvec_folder=$out_dir/empl_corrected_bvec
emp_bval_folder=$out_dir/emp_corrected_bval
mkdir $emp_bvec_folder && $emp_bval_folder
mv $out_dir/bvec_00* $emp_bvec_folder && mv $out_dir/bval_00* $emp_bval_folder

# Compute FA/MD/PEV of empirical correction
out_name=emp_corrected
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "dti_voxel_fit(uncorrected_sig,bvec_folder,bval_folder,mask_path, out_dir, out_name)"
com

