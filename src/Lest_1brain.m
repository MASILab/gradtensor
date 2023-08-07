
img = double(niftiread('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs001/ses-s1Bx2/dwi/sub-cIVs001_ses-s1Bx2_acq-b2000n56r21x21x22peAPP_run-111_dwi.nii'));
dt = niftiread('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs001/ses-s1Bx2/dwi_metric/dt_b2000_run-1.nii.gz');
g = importdata('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs001/ses-s1Bx2/dwi/sub-cIVs001_ses-s1Bx2_acq-b2000n56r21x21x22peAPP_run-111_dwi.bvec');
b = importdata('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs001/ses-s1Bx2/dwi/sub-cIVs001_ses-s1Bx2_acq-b2000n56r21x21x22peAPP_run-111_dwi.bval');
b0 = niftiread('/nfs/masi/kanakap/projects/estimate/MASiVar_kids/sub-cIVs001/ses-s1Bx2/dwi_metric/meanb0_b2000run-1.nii.gz');
%%
i = 55; j =73; k = 32;
L_pred = test_Lvol_wholeimg1(dt(i,j,k,:),g,b,img(i,j,k,:),b0(i,j,k));


%%

%L_pred(i,j,k,:) = test_Lvol_wholeimg2(DT_back_to_bore_stack(i,j,k,:,:,:),bvec_all,bval_all_nonb0,dwi_all(i,j,k,:,:),b0_all(i,j,k,:));

%%

% true_L = '/nfs/masi/kanakap/projects/estimate/L.nii';
% Vref = spm_vol(true_L);
