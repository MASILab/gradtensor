# lr pipeline

# Inputs 1 -> sub/sess, 2 -> dwi, 3 -> bvec, 4 -> bval, 5 -> SNR, 6 -> direction, 7 -> run #
echo "STARTING PIPELINE";
echo "$1 $2 $3 $4 $5 $6 $7 $8"
mkdir /tmp/decompress;
dwi_path_gz=/nobackup/p_masi_brain_map/kanakap/LR/MASIver_kids_prequal-v1.0.0/$1/dwi/$2.gz;
gzip -dc $dwi_path_gz > /tmp/decompress/$2;
dwi_path=/tmp/decompress/$2;
org_bvec_path=/nobackup/p_masi_brain_map/kanakap/LR/MASIver_kids_prequal-v1.0.0/$1/dwi/$4;
org_bval_path=/nobackup/p_masi_brain_map/kanakap/LR/MASIver_kids_prequal-v1.0.0/$1/dwi/$3;
#mask_path=/nfs/masi/kanakap/projects/LR/masivar_input/$4
initial_SNR=inf;
out_dir=/nobackup/p_masi_brain_map/kanakap/LR/MASiVar_kids/$1/tracto_ip_$5_$6;
mkdir $out_dir;
uncorr_out_name=uncorrected;
Nest_uncorr_out_name=uncorrected_Nest;
Lest_uncorr_out_name=uncorrected_Lest;

uncorrected_sig=$out_dir/uncorrected_NLest_sig.nii;
Nuncorrected_sig=$out_dir/uncorrected_Nest_sig.nii;
Luncorrected_sig=$out_dir/uncorrected_Lest_sig.nii;

Limg_file=/scratch/kanakap/LR/L.nii.gz;

# Reconstruct corrput signal
echo "RUNNING LR corrpution";
matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages');addpath(genpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/pkgs/spm_read_nii'));addpath(genpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/pkgs/spm_reslice'));addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/pkgs/sh_basis');addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/pkgs/niftilib/');addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/pkgs/revised_matlab_functions/');disp('Reconstruct corrput signal');/nobackup/p_masi_brain_map/kanakap/LR/scripts/compute_noise_corrput_signal_v1('$dwi_path','$mask_path','$out_dir','$uncorr_out_name','$Limg_file','$org_bvec_path','$org_bval_path',$initial_SNR);disp('Compute FA/MD/PEV of LR+noise signal'); exit";

# Make the input and output dir for Francois speical
rm -r /tmp/decompress;
cp $org_bvec_path $out_dir;
cp $org_bval_path $out_dir;
anat_dir=/nobackup/p_masi_brain_map/kanakap/LR/MASiVar_kids/$1/anat;
seg_dir=/nobackup/p_masi_brain_map/kanakap/LR/MASiVar_kids/$1/slant_output/FinalResult;
anat_file=$(ls $anat_dir);
seg_file=$(ls $seg_dir);
cp $anat_dir/* $out_dir;
cp $seg_dir/* $out_dir;
francois_out=/nobackup/p_masi_brain_map/kanakap/LR/MASiVar_kids/$1/tracto_op_$5_$6_est;
mkdir $francois_out;
gzip -f $out_dir/uncorrected_est_sig.nii;

# Francois speical command
JOBDIR=/tmp/francois_$5_$6_est;
mkdir $JOBDIR
INDIR=$out_dir;
OUTDIR=$francois_out;

echo "RUNNING Francois speical est";
echo "singularity run --home $JOBDIR --bind $JOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $OUTDIR:/OUTPUTS --bind $JOBDIR:/TMP /nobackup/p_masi_brain_map/kanakap/LR/singularity_francois_special_v1.sif $7 $8 $anat_file uncorrected_est_sig.nii.gz $3 $4 8 "0 $5" "0 $5" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52"";
singularity run --home $JOBDIR --bind $JOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $OUTDIR:/OUTPUTS --bind $JOBDIR:/TMP /scratch/kanakap/LR/singularity_francois_special_v1.sif $7 $8 $anat_file uncorrected_est_sig.nii.gz $3 $4 8 "0 $5" "0 $5" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52";

echo "RUNNING Francois speical Lest";
gzip -f $out_dir/uncorrected_Lest_sig.nii;
Lfrancois_out=/nobackup/p_masi_brain_map/kanakap/LR/MASiVar_kids/$1/tracto_op_$5_$6_Lest;
mkdir $Lfrancois_out;
LJOBDIR=/tmp/francois_$5_$6_Lest
mkdir $LJOBDIR;
LOUTDIR=$Lfrancois_out;
echo "singularity run --home $LJOBDIR --bind $LJOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $LOUTDIR:/OUTPUTS --bind $LJOBDIR:/TMP /nobackup/p_masi_brain_map/kanakap/LR/singularity_francois_special_v1.sif $7 $8 $anat_file uncorrected_Lest_sig.nii.gz $3 $4 8 "0 $5" "0 $5" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52"";
singularity run --home $LJOBDIR --bind $LJOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $LOUTDIR:/OUTPUTS --bind $LJOBDIR:/TMP /scratch/kanakap/LR/singularity_francois_special_v1.sif $7 $8 $anat_file uncorrected_Lest_sig.nii.gz $3 $4 8 "0 $5" "0 $5" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52";
