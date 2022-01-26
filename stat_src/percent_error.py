import nibabel as nib
import numpy as np
import scipy.io as sio

true_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_fa.nii').get_fdata()
corpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_noflip/Lest_fa.nii').get_fdata()
sm_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_noflip/Lest_corr_sm_fa.nii').get_fdata()
bx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_noflip/Lest_corr_bx_fa.nii').get_fdata()

true_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii').get_fdata()
corpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_noflip/Lest_md.nii').get_fdata()
sm_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_noflip/Lest_corr_sm_md.nii').get_fdata()
bx_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_noflip/Lest_corr_bx_md.nii').get_fdata()

true_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_primary_eigvec.nii').get_fdata()
corpt_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/OUTPUTS_estimates_noflip/Lest_primary_eigvec.nii').get_fdata()
sm_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_sm_noflip/Lest_corr_sm_primary_eigvec.nii').get_fdata()
bx_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posA/ISMRM_bx_noflip/Lest_corr_bx_primary_eigvec.nii').get_fdata()

btrue_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_fa.nii').get_fdata()
bcorpt_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_estimates_noflip/Lest_fa.nii').get_fdata()
bsm_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_sm_noflip/Lest_corr_sm_fa.nii').get_fdata()
bbx_corr_fa = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_bx_noflip/Lest_corr_bx_fa.nii').get_fdata()

btrue_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_md.nii').get_fdata()
bcorpt_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_estimates_noflip/Lest_md.nii').get_fdata()
bsm_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_sm_noflip/Lest_corr_sm_md.nii').get_fdata()
bbx_corr_md = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_bx_noflip/Lest_corr_bx_md.nii').get_fdata()

btrue_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_future_fieldmap/p_3tb_posB_mask_primary_eigvec.nii').get_fdata()
bcorpt_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/OUTPUTS_estimates_noflip/Lest_primary_eigvec.nii').get_fdata()
bsm_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_sm_noflip/Lest_corr_sm_primary_eigvec.nii').get_fdata()
bbx_corr_pev = nib.load('/home/local/VANDERBILT/kanakap/gradtensor_data/10_29_2019_human_repositioned/3tb/posB/ISMRM_bx_noflip/Lest_corr_bx_primary_eigvec.nii').get_fdata()

def rmse_mat(imga, imgb):
    err = np.nanmean((imga - imgb) ** 2)
          #den = float(96 * 96 * 68)
    err = np.sqrt(err)
    return err

def rmse_voxelwise(imga, imgb):
    rmse_img = np.zeros((96,96,68))
    for x in range(96):
        for y in range(96):
            for z in range(68):
                err = np.mean((imga[x,y,z] - imgb[x,y,z]) ** 2)
                err = np.sqrt(err)
                rmse_img[x,y,z] = err 
    return rmse_img

def angular_error(PEa, PEb, halfPi=True):
    # PEa = preprocessing.normalize(PEa, norm='l1', axis=1)
    # PEb = preprocessing.normalize(PEb, norm='l1', axis=1)
    chord = np.square(PEa[..., 0] - PEb[..., 0]) + \
    np.square(PEa[..., 1] - PEb[..., 1]) + \
    np.square(PEa[..., 2] - PEb[..., 2])
    chord = np.sqrt(chord)
    ang = 2 * np.real(np.arcsin(chord/2))
    if halfPi:
        ang[ang > (np.pi/2)] = np.pi - ang[ang > (np.pi/2)]
    return np.degrees(ang)

#bins = np.histogram(np.hstack((np.nan_to_num(after_sm_fa),np.nan_to_num(before))),bins=20)[1]
#plt.hist(after_sm_fa.ravel(),bins=bins,color = ['r'], label = 'after simple correction', alpha = 0.5)
#plt.hist(before.ravel(),bins=bins,color = ['b'], label = 'before correction', alpha = 0.5)
#plt.legend()
#plt.ylabel('number of voxels')
#plt.xlabel('rmse fa')
#plt.title('Histogram of FA before and simple correction')
#plt.show()

#true_lower_after = []
#for i in range(96):
#    for j in range(96):
#        for k in range(68):
#            if before[i,j,k] > after_sm_fa[i,j,k]:
#                true_lower_after.append(true_fa[i,j,k])
#plt.scatter(range(len(LR_higher)),LR_higher)

#LR_higher = []
#for i in range(96):
#    for j in range(96):
#        for k in range(68):
#            if before[i,j,k] > after_sm_fa[i,j,k]:
#                L_mat = np.squeeze(vL[:,:,x,y,z])
#                LR_higher.append((np.linalg.det(L_mat))

percent_change_fa_sm = 1 - (rmse_mat(true_fa, sm_corr_fa) / rmse_mat(true_fa, corpt_fa))
percent_change_fa_bx = 1 - (rmse_mat(true_fa, bx_corr_fa) / rmse_mat(true_fa, corpt_fa))
print('fa full', percent_change_fa_bx)
print('fa simple', percent_change_fa_sm)
percent_change_md_sm = 1 - (rmse_mat(true_md, sm_corr_md) / rmse_mat(true_md, corpt_md))
percent_change_md_bx = 1 - (rmse_mat(true_md, bx_corr_md) / rmse_mat(true_md, corpt_md))
print('md full', percent_change_md_bx)
print('md simple', percent_change_md_sm)


bpercent_change_fa_sm = 1 - (rmse_mat(btrue_fa, bsm_corr_fa) / rmse_mat(btrue_fa, bcorpt_fa))
bpercent_change_fa_bx = 1 - (rmse_mat(btrue_fa, bbx_corr_fa) / rmse_mat(btrue_fa, bcorpt_fa))
print('fa full pos b', bpercent_change_fa_bx)
print('fa simple pos b', bpercent_change_fa_sm)
bpercent_change_md_sm = 1 - (rmse_mat(btrue_md, bsm_corr_md) / rmse_mat(btrue_md, bcorpt_md))
bpercent_change_md_bx = 1 - (rmse_mat(btrue_md, bbx_corr_md) / rmse_mat(btrue_md, bcorpt_md))
print('md full pos b', bpercent_change_md_bx)
print('md simple pos b', bpercent_change_md_sm)

percent_change_pev_sm = 1 - (angular_error(true_pev, sm_corr_pev) / angular_error(true_pev, corpt_pev))
percent_change_pev_bx = 1 - (angular_error(true_pev, bx_corr_pev) / angular_error(true_pev, corpt_pev))
print('pev full ', np.nanmean(percent_change_pev_bx))
print('pev simple ', np.nanmean(percent_change_pev_sm))

bpercent_change_pev_sm = 1 - (angular_error(btrue_pev, bsm_corr_pev) / angular_error(btrue_pev, bcorpt_pev))
bpercent_change_pev_bx = 1 - (angular_error(btrue_pev, bbx_corr_pev) / angular_error(btrue_pev, bcorpt_pev))
print('md full pos b', np.nanmean(bpercent_change_pev_bx))
print('md simple pos b', np.nanmean(bpercent_change_pev_sm))

perc_err = 100* ((true_fa - corpt_fa) / true_fa)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)

perc_err = 100* ((true_md - corpt_md) / true_md)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)

perc_err = 100* ((true_fa - bx_corr_fa) / true_fa)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)

perc_err = 100* ((true_md - bx_corr_md) / true_md)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)

perc_err = 100* ((true_fa - sm_corr_fa) / true_fa)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)

perc_err = 100* ((true_md - sm_corr_md) / true_md)
abs_perc_err = np.abs(perc_err)
np.nanmean(abs_perc_err)
