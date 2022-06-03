# lr pipeline

# Inputs 1 -> sub/sess, 2 -> dwi, 3 -> bvec, 4 -> bval, 5 -> SNR, 6 -> direction, 7 -> run #
mkdir /tmp/decompress
dwi_path_gz=/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/$1/dwi/$2.gz
gzip -dc $dwi_path_gz > /tmp/decompress/$2
dwi_path=/tmp/decompress/$2
org_bvec_path=/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/$1/dwi/$3
org_bval_path=/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/$1/dwi/$4
#mask_path=/nfs/masi/kanakap/projects/LR/masivar_input/$4
initial_SNR=$5
out_dir=/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/$1/tracto_ip_$6_$7

mkdir $out_dir
uncorr_out_name=uncorrected
Nest_uncorr_out_name=uncorrected_Nest
Lest_uncorr_out_name=uncorrected_Lest

uncorrected_sig=$out_dir/uncorrected_NLest_sig.nii
Nuncorrected_sig=$out_dir/uncorrected_Nest_sig.nii
Luncorrected_sig=$out_dir/uncorrected_Lest_sig.nii

Limg_file=/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz

# Reconstruct corrput signal
~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_read_nii'));addpath(genpath('/home/local/VANDERBILT/kanakap/gradtensor/external/spm_reslice'));addpath('/home/local/VANDERBILT/kanakap/gradtensor/src/sh_basis');addpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/');addpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/');disp('Reconstruct corrput signal');compute_noise_corrput_signal_v1('$dwi_path','$mask_path','$out_dir','$uncorr_out_name','$Limg_file','$org_bvec_path','$org_bval_path',$initial_SNR);disp('Compute FA/MD/PEV of LR+noise signal'); exit"

rm -r /tmp/decompress
cp $org_bvec_path $out_dir
cp $org_bval_path $out_dir
ln -s /nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/$1/anat/* $out_dir
