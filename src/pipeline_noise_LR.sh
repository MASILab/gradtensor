# lr pipeline

# Inputs
#dwi_path=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/INPUTS/dwmri.nii
dwi_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.nii
bvec_folder=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bvec/
bval_folder=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/corrected_bval/
org_bvec_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bvec
org_bval_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/sub-cIIs00_ses-s1Bx1_acq-b1000n32r25x25x25peAPP_run-112_dwi.bval
mask_path=/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii

out_dir=/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1
#rL_path=/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/L_resamp.nii

mkdir $out_dir
uncorr_out_name=uncorrected
Nest_uncorr_out_name=uncorrected_Nest
Lest_uncorr_out_name=uncorrected_Lest

initial_SNR=10
uncorrected_sig=$out_dir/uncorrected_NLest_sig.nii
Nuncorrected_sig=$out_dir/uncorrected_Nest_sig.nii
Luncorrected_sig=$out_dir/uncorrected_Lest_sig.nii

Limg_file=/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz
#refimg_file=$uncorrected_sig
approx_out_name=approx_corrected
Lapprox_out_name=approx_Lcorrected
Napprox_out_name=approx_Ncorrected

approx_corrected=$out_dir/approx_corrected_sig.nii
Lapprox_corrected=$out_dir/approx_Lcorrected_sig.nii
Napprox_corrected=$out_dir/approx_Ncorrected_sig.nii

emp_out_dir=$out_dir/emp
Nemp_out_dir=$out_dir/Nemp
Lemp_out_dir=$out_dir/Lemp
mkdir $emp_out_dir
mkdir $Nemp_out_dir
mkdir $Lemp_out_dir

# Reconstruct corrput signal
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_read_nii'));addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_reslice'));addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/sh_basis');addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');disp('Reconstruct corrput signal');compute_noise_corrput_signal('$dwi_path','$bvec_folder','$bval_folder','$mask_path','$out_dir','$uncorr_out_name','$Limg_file','$org_bvec_path','$org_bval_path',$initial_SNR);disp('Compute FA/MD/PEV of LR+noise signal');dti_voxel_fit_sig('$uncorrected_sig','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$uncorr_out_name');disp('Compute FA/MD/PEV of noise signal');dti_voxel_fit_sig('$Nuncorrected_sig','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$Nest_uncorr_out_name');disp('Compute FA/MD/PEV of LR signal');dti_voxel_fit_sig('$Luncorrected_sig','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$Lest_uncorr_out_name');disp('Approximate correction LR+noise');compute_approx_corr_dwi('$Limg_file','$uncorrected_sig','$org_bval_path','$org_bvec_path','$mask_path','$approx_out_name','$out_dir');disp('Approximate correction noise');compute_approx_corr_dwi('$Limg_file','$Nuncorrected_sig','$org_bval_path','$org_bvec_path','$mask_path','$Napprox_out_name','$out_dir'); disp('Approximate correction LR');compute_approx_corr_dwi('$Limg_file','$Luncorrected_sig','$org_bval_path','$org_bvec_path','$mask_path','$Lapprox_out_name','$out_dir'); disp('Compute FA/MD/PEV of LR+noise approximate correction');dti_voxel_fit_sig('$approx_corrected','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$approx_out_name');disp('Compute FA/MD/PEV of noise approximate correction');dti_voxel_fit_sig('$Napprox_corrected','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$Napprox_out_name');disp('Compute FA/MD/PEV of LR approximate correction');dti_voxel_fit_sig('$Lapprox_corrected','$org_bvec_path','$org_bval_path','$mask_path','$out_dir','$Lapprox_out_name');disp('Empirical correction LR+noise');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$uncorrected_sig','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$emp_out_dir');disp('Empirical correction noise');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$Nuncorrected_sig','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$Nemp_out_dir');disp('Empirical correction LR');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$Luncorrected_sig','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$Lemp_out_dir'); exit"

emp_bvec_folder=$emp_out_dir/emp_corrected_bvec
emp_bval_folder=$emp_out_dir/emp_corrected_bval
mkdir $emp_bvec_folder
mkdir $emp_bval_folder
mv $emp_out_dir/bvec_00* $emp_bvec_folder
mv $emp_out_dir/bval_00* $emp_bval_folder

Nemp_bvec_folder=$Nemp_out_dir/emp_Ncorrected_bvec
Nemp_bval_folder=$Nemp_out_dir/emp_Ncorrected_bval
mkdir $Nemp_bvec_folder
mkdir $Nemp_bval_folder
mv $Nemp_out_dir/bvec_00* $Nemp_bvec_folder
mv $Nemp_out_dir/bval_00* $Nemp_bval_folder

Lemp_bvec_folder=$Lemp_out_dir/emp_Lcorrected_bvec
Lemp_bval_folder=$Lemp_out_dir/emp_Lcorrected_bval
mkdir $Lemp_bvec_folder
mkdir $Lemp_bval_folder
mv $Lemp_out_dir/bvec_00* $Lemp_bvec_folder
mv $Lemp_out_dir/bval_00* $Lemp_bval_folder

emp_out_name=emp_corrected
Nemp_out_name=Nemp_corrected
Lemp_out_name=Lemp_corrected

~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('../external/spm_read_nii'));addpath(genpath('../external/spm_reslice'));addpath('sh_basis');addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');disp('Compute FA/MD/PEV of LR+noise empirical correction');dti_voxel_fit_masiver('$uncorrected_sig','$emp_bvec_folder','$emp_bval_folder','$mask_path','$emp_out_dir','$emp_out_name','$org_bval_path','$org_bvec_path');disp('Compute FA/MD/PEV of noise empirical correction');dti_voxel_fit_masiver('$Nuncorrected_sig','$Nemp_bvec_folder','$Nemp_bval_folder','$mask_path','$Nemp_out_dir','$Nemp_out_name','$org_bval_path','$org_bvec_path');disp('Compute FA/MD/PEV of LR empirical correction');dti_voxel_fit_masiver('$Luncorrected_sig','$Lemp_bvec_folder','$Lemp_bval_folder','$mask_path','$Lemp_out_dir','$Lemp_out_name','$org_bval_path','$org_bvec_path');exit"

