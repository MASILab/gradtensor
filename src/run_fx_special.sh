

. final_pipeline_LR_prequal.sh $1 $2.nii $3 $4 $5 $6


JOBDIR=/tmp/francois_run1; INDIR=/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_ip_1000_113;OUTDIR=/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/sub-cIVs001/ses-s1Bx2/tracto_op_1000_113;singularity run --home $JOBDIR --bind $JOBDIR:/tmp --containall --cleanenv --bind $INDIR:/INPUTS --bind $OUTDIR:/OUTPUTS --bind $JOBDIR:/TMP /nfs/masi/kanakap/singularity_francois_special_v1.sif sub-cIVs001 ses-s1Bx2 sub-cIVs001_ses-s1Bx2_acq-r10x10x10_T1w.nii.gz uncorrected_est_sig.nii.gz sub-cIVs001_ses-s1Bx2_acq-b1000n40r21x21x22peAPP_run-113_dwi.bval sub-cIVs001_ses-s1Bx2_acq-b1000n40r21x21x22peAPP_run-113_dwi.bvec 8 "0 1000" "0 1000" 1 wm 1 wm prob 27 0.4 20 sub-cIVs001_ses-s1Bx2_acq-r10x10x10_T1w_seg.nii.gz "4 40 41 44 45 51 52"
