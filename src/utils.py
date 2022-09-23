import torch
import logging
import warnings
import numpy as np
from dipy.core.sphere import Sphere
from dipy.reconst.shm import sf_to_sh, sh_to_sf
from dipy.core.gradients import gradient_table_from_bvals_bvecs
from scilpy.reconst.raw_signal import compute_sh_coefficients
from scilpy.utils.bvec_bval_tools import (check_b0_threshold, identify_shells,
                                          is_normalized_bvecs, normalize_bvecs)



def compute_dwi_attenuation(dwi_weights: np.ndarray, b0: np.ndarray):
    """ Compute signal attenuation by dividing the dwi signal with the b0.
    Parameters:
    -----------
    dwi_weights : np.ndarray of shape (X, Y, Z, #gradients)
        Diffusion weighted images.
    b0 : np.ndarray of shape (X, Y, Z)
        B0 image.
    Returns
    -------
    dwi_attenuation : np.ndarray
        Signal attenuation (Diffusion weights normalized by the B0).
    """
    b0 = b0[..., None]  # Easier to work if it is a 4D array.

    # Make sure that, in every voxels, weights are lower in the b0. Should
    # always be the case, but with the noise we never know!
    erroneous_voxels = np.any(dwi_weights > b0, axis=3)
    nb_erroneous_voxels = np.sum(erroneous_voxels)
    if nb_erroneous_voxels != 0:
        logging.info("# of voxels where `dwi_signal > b0` in any direction: "
                     "{}".format(nb_erroneous_voxels))
        dwi_weights = np.minimum(dwi_weights, b0)

    # Compute attenuation
    dwi_attenuation = dwi_weights / b0

    # Make sure we didn't divide by 0.
    dwi_attenuation[np.logical_not(np.isfinite(dwi_attenuation))] = 0.

    return dwi_attenuation


# DWI TO SF (With the voxelwise bvec and bval)
def reconstruct_signal_at_voxel(i,j,k,n,og_bvec,og_bval,bvec_stack,bval_stack,dwi_hat,sh_order):
        warnings.filterwarnings("ignore")
        dwi = n[i][j][k]
        og_gradient_table = gradient_table_from_bvals_bvecs(og_bval, og_bvec)
        vec = bvec_stack[i,j,k,:,:]
        val = bval_stack[i,j,k,:]
        gradient_table = gradient_table_from_bvals_bvecs(val, vec)
        sh_order=sh_order
        basis_type='tournier07'
        smooth=0.00
        use_attenuation=True
        force_b0_threshold=False
        mask=None
        sphere=None

        # Extracting infos
        b0_mask = gradient_table.b0s_mask
        bvecs = gradient_table.bvecs
        bvals = gradient_table.bvals

        dwi = np.reshape(dwi,[1,1,1,bvals.shape[0]])

        if not is_normalized_bvecs(bvecs):
                logging.warning("Your b-vectors do not seem normalized...")
                bvecs = normalize_bvecs(bvecs)

        b0_threshold = check_b0_threshold(force_b0_threshold, bvals.min())

        # Ensure that this is on a single shell.
        shell_values, _ = identify_shells(bvals)
        shell_values.sort()
        # if shell_values.shape[0] != 2 or shell_values[0] > b0_threshold:
        #     raise ValueError("Can only work on single shell signals.")

        # Keeping b0-based infos
        bvecs = bvecs[np.logical_not(b0_mask)]
        weights = dwi[..., np.logical_not(b0_mask)]

        # scale singal with bval correction
        b0 = dwi[..., b0_mask].mean(axis=3)
        norm_gg = np.divide(bvals[np.logical_not(b0_mask)] , og_bval[np.logical_not(b0_mask)])
        weights_scaled = b0 * np.exp (np.divide( (np.log (np.divide(weights,b0)) ) , norm_gg))

        # Compute attenuation using the b0.
        if use_attenuation:
                weights_scaled = compute_dwi_attenuation(weights_scaled, b0)

        # # Get cartesian coords from bvecs
        sphere = Sphere(xyz=bvecs)

        # SF TO SH
        # Fit SH
        sh = sf_to_sh(weights_scaled, sphere, sh_order, basis_type, smooth=smooth)

        # Apply mask
        if mask is not None:
                sh *= mask[..., None]

        # Reconstructing DWI
        # SH to SF
        og_bvecs = og_gradient_table.bvecs

        if not is_normalized_bvecs(og_bvecs):
                logging.warning("Your b-vectors do not seem normalized...")
                og_bvecs = normalize_bvecs(og_bvecs)

        og_bvecs = og_bvecs[np.logical_not(b0_mask)]

        og_sphere = Sphere(xyz=og_bvecs)

        sf = sh_to_sf(sh, og_sphere, sh_order=sh_order, basis_type=basis_type)

        # SF TO DWI (inverse of compute_dwi_attenuation) here weights_hat is DWI with bvec corrected
        b0 = b0[..., None]
        weights_hat = sf * b0
        dwi_hat[i,j,k,:] = weights_hat

def val_pk(data,og_bval,og_bvec,sh_order):
    og_gradient_table = gradient_table_from_bvals_bvecs(og_bval, og_bvec)
    pk_sh = compute_sh_coefficients(data,og_gradient_table,sh_order=sh_order,basis_type='tournier07',use_attenuation=True,smooth=0.00)
    return pk_sh

def val_emp(i,j,k,n,og_bval,bvec_stack,bval_stack,emp_sh,sh_order):
        warnings.filterwarnings("ignore")
        dwi = n[i][j][k]
        vec = bvec_stack[i,j,k,:,:]
        val = bval_stack[i,j,k,:]
        gradient_table = gradient_table_from_bvals_bvecs(val, vec)
        sh_order=sh_order
        basis_type='tournier07'
        smooth=0.00
        use_attenuation=True
        force_b0_threshold=False
        mask=None
        sphere=None

        # Extracting infos
        b0_mask = gradient_table.b0s_mask
        bvecs = gradient_table.bvecs
        bvals = gradient_table.bvals
        
        dwi = np.reshape(dwi,[1,1,1,bvals.shape[0]])

        if not is_normalized_bvecs(bvecs):
                logging.warning("Your b-vectors do not seem normalized...")
                bvecs = normalize_bvecs(bvecs)

        b0_threshold = check_b0_threshold(force_b0_threshold, bvals.min())

        # Ensure that this is on a single shell.
        shell_values, _ = identify_shells(bvals)
        shell_values.sort()
        # if shell_values.shape[0] != 2 or shell_values[0] > b0_threshold:
        #     raise ValueError("Can only work on single shell signals.")

        # Keeping b0-based infos
        bvecs = bvecs[np.logical_not(b0_mask)]
        weights = dwi[..., np.logical_not(b0_mask)]

        b0 = dwi[..., b0_mask].mean(axis=3)
        norm_gg = np.divide(bvals[np.logical_not(b0_mask)] , og_bval[np.logical_not(b0_mask)])
        weights_scaled = b0 * np.exp (np.divide( (np.log (np.divide(weights,b0)) ) , norm_gg))

        # Compute attenuation using the b0.
        if use_attenuation:
                weights_scaled = compute_dwi_attenuation(weights_scaled, b0)

        # # Get cartesian coords from bvecs # from here cut debugging
        sphere = Sphere(xyz=bvecs)

        # SF TO SH
        # Fit SH
        sh = sf_to_sh(weights_scaled, sphere, sh_order, basis_type, smooth=smooth)
        emp_sh[i,j,k,:] = sh


def angularCorrCoeff(p_v, q_v):

    p = p_v.shape[-1]
    q = q_v.shape[-1]
    nCoeffs = min([p, q])

    p_v = p_v[..., 0:nCoeffs]
    q_v = q_v[..., 0:nCoeffs]

    zmp_v = p_v[..., 1:nCoeffs]
    zmq_v = q_v[..., 1:nCoeffs]

    np_v = np.rollaxis(zmp_v, -1)/ (np.sqrt(np.sum(zmp_v * np.conj(zmp_v), axis=-1)))
    nq_v = np.rollaxis(zmq_v, -1)/ (np.sqrt(np.sum(zmq_v * np.conj(zmq_v), axis=-1)))
    
    np_v = np.rollaxis(np_v,0,len(np_v.shape))
    nq_v = np.rollaxis(nq_v,0,len(nq_v.shape))

    acc = np.sum(np_v * np.conj(nq_v), axis=-1)

    return acc