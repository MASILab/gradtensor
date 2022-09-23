# lr pipeline

# Inputs 1 -> sub/sess, 2 -> dwi, 3 -> bval, 4 -> bvec, 5 -> run, 6 -> subj, 7 -> sess, 8 -> root dir
echo "STARTING PIPELINE ON PREPROCESSED DATA (next lr correction) ";
echo "$1 $2 $3 $4 $5 $6 $7 $8"
mkdir /tmp/decompress;
dwi_path_gz=$8/$1/prequal_dwi_cat/$2.gz;
gzip -dc $dwi_path_gz > /tmp/decompress/$2;
dwi_path=/tmp/decompress/$2;
org_bval_path=$8/$1/prequal_dwi_cat/$3;
org_bvec_path=$8/$1/prequal_dwi_cat/$4;
out_dir=$8/$1/tracto_ip_lr_corr_$5;
mkdir $out_dir;

Luncorrected_sig=$dwi_path; # should be dwi from prequal - this is uncorrected
Lemp_out_dir=$out_dir/Lemp
mkdir $Lemp_out_dir

#Limg_file=$6/$1/L.nii;
Limg_file=/home/local/VANDERBILT/kanakap/gradtensor_data/fieldmaps/3tb_future_fieldmap/OUTPUTS/L.nii.gz;

# Empircal correction
#pkecho "RUNNING LR EMP CORRECTION";
#matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages'),addpath(genpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts')),addpath(genpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/spm_read_nii')),addpath(genpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/spm_reslice')),addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/sh_basis'),addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/niftilib/'),addpath('/nobackup/p_masi_brain_map/kanakap/LR/scripts/revised_matlab_functions/'),disp('Empirical correction LR');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$Luncorrected_sig','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$Lemp_out_dir'); exit"

#pk~/MATLAB_2017a_install/bin/matlab -nodisplay -nosplash -nodesktop -r "disp('adding packages'),addpath(genpath('/nfs/masi/kanakap/xnat_apps/masimatab/trunk/xnatspiders/matlab/justinlib_v1_7_0/niftilib/')),addpath(genpath('/home/local/VANDERBILT/kanakap/XNAT/TemporalLobe/revised_matlab_functions/')),addpath(genpath('../external/spm_read_nii')),addpath('../external/spm_reslice'),disp('Empirical correction LR');apply_gradtensor_to_b('Limg_file','$Limg_file','refimg_file','$Luncorrected_sig','bval_file','$org_bval_path','bvec_file','$org_bvec_path','out_dir','$Lemp_out_dir'); exit"

# Move to the correct folders
Lemp_out_dir=$out_dir/Lemp;
mkdir $Lemp_out_dir;
Lemp_bvec_folder=$Lemp_out_dir/emp_Lcorrected_bvec;
Lemp_bval_folder=$Lemp_out_dir/emp_Lcorrected_bval;
mkdir $Lemp_bvec_folder;
mkdir $Lemp_bval_folder;
mv $Lemp_out_dir/bvec_* $Lemp_bvec_folder;
mv $Lemp_out_dir/bval_* $Lemp_bval_folder;
Lemp_out_name=emp;
emp_sig=$out_dir/emp_Lcorrected_sig.nii.gz;
echo $emp_sig
# Voxelwise to signal 
echo "RUNNING LR EMP CORRECTION SIGNAL ESTIMATION";
#source /nobackup/p_masi_brain_map/kanakap/LR/scripts/py38/bin/activate
#pkpython lr_sh_to_sig.py $dwi_path $org_bvec_path $org_bval_path $Lemp_bvec_folder $Lemp_bval_folder $out_dir $Lemp_out_name
python lr_sh_val.py $dwi_path $org_bvec_path $org_bval_path $Lemp_bvec_folder $Lemp_bval_folder $emp_sig

<<com
# Make the input and output dir for Francois speical
cp $org_bvec_path $out_dir;
cp $org_bval_path $out_dir;
anat_dir=$8/$1/anat;
seg_dir=$8/$1/slant_output/FinalResult;
anat_file=$(ls $anat_dir);
seg_file=$(ls $seg_dir);
cp $anat_dir/* $out_dir;
cp $seg_dir/* $out_dir;
francois_out=$8/$1/tracto_op_lr_emp_corr_$5;
mkdir $francois_out;
#gzip -f $out_dir/emp_Lcorrected_sig.nii;

# Francois speical command
JOBDIR=/tmp/francois_$6_$7_$5_corr;
mkdir $JOBDIR;
INDIR=$out_dir;
OUTDIR=$francois_out;

echo "RUNNING Francois speical on corrected Lr";
echo "singularity run --home $JOBDIR --bind $JOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $OUTDIR:/OUTPUTS --bind $JOBDIR:/TMP /nfs/masi/kanakap/singularity_francois_special_v1.sif $6 $7 $anat_file emp_Lcorrected_sig.nii.gz $3 $4 8 "0 1000" "0 2000" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52"";

singularity run --home $JOBDIR --bind $JOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $OUTDIR:/OUTPUTS --bind $JOBDIR:/TMP /nfs/masi/kanakap/singularity_francois_special_v1.sif $6 $7 $anat_file emp_Lcorrected_sig.nii.gz $3 $4 8 "0 1000" "0 2000" 5 wm 5 wm prob 27 0.4 20 $seg_file "4 40 41 44 45 51 52";

rm -r $JOBDIR
com
