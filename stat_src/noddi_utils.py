import nibabel as nib 
import numpy as np

def get_data(scalar_file):
    root = '/home/local/VANDERBILT/kanakap/gradtensor/src/noddi'
    mask_root = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/'
    bammer1 = nib.load(root + '/noddi_output_corr1/noddi_' + scalar_file).get_fdata()
    bammer2 = nib.load(root + '/noddi_output_corr2/noddi_' + scalar_file).get_fdata()
    bammer3 = nib.load(root + '/noddi_output_corr3/noddi_' + scalar_file).get_fdata()
    bammer4 = nib.load(root + '/noddi_output_corr4/noddi_' + scalar_file).get_fdata()
    bammer5 = nib.load(root + '/noddi_output_corr5/noddi_' + scalar_file).get_fdata()

    prop1 = nib.load(root + '/noddi_output_proposed1/noddi_' + scalar_file).get_fdata()
    prop2 = nib.load(root + '/noddi_output_proposed2/noddi_' + scalar_file).get_fdata()
    prop3 = nib.load(root + '/noddi_output_proposed3/noddi_' + scalar_file).get_fdata()
    prop4 = nib.load(root + '/noddi_output_proposed4/noddi_' + scalar_file).get_fdata()
    prop5 = nib.load(root + '/noddi_output_proposed5/noddi_' + scalar_file).get_fdata()

    bammer1[bammer1 == 0] = 'nan'
    bammer2[bammer2 == 0] = 'nan'
    bammer3[bammer3 == 0] = 'nan'
    bammer4[bammer4 == 0] = 'nan'
    bammer5[bammer5 == 0] = 'nan'
    prop1[prop1 == 0] = 'nan'
    prop2[prop2 == 0] = 'nan'
    prop3[prop3 == 0] = 'nan'
    prop4[prop4 == 0] = 'nan'
    prop5[prop5 == 0] = 'nan'
    
    mask1 = nib.load(mask_root + 'sub-cIVs001/ses-s1Bx2/fast_reg1/LABELS_IN_DWI.nii.gz').get_fdata()
    mask2 = nib.load(mask_root + 'sub-cIVs002/ses-s1Bx2/fast_reg1/LABELS_IN_DWI.nii.gz').get_fdata()
    mask3 = nib.load(mask_root + 'sub-cIVs006/ses-s1Bx2/fast_reg1/LABELS_IN_DWI.nii.gz').get_fdata()
    mask4 = nib.load(mask_root + 'sub-cIVs007/ses-s1Bx2/fast_reg1/LABELS_IN_DWI.nii.gz').get_fdata()
    mask5 = nib.load(mask_root + 'sub-cIVs005/ses-s1Bx3/fast_reg1/LABELS_IN_DWI.nii.gz').get_fdata()

    if scalar_file != 'fiso.nii':
        bammer1 = np.where(mask1 != 1.0,bammer1, np.nan) # wm gm 
        bammer1 = np.where(mask1 != 1.0, bammer1,np.nan)
        bammer2 = np.where(mask2 != 1.0, bammer2,np.nan)
        bammer3 = np.where(mask3 != 1.0, bammer3,np.nan)
        bammer4 = np.where(mask4 != 1.0, bammer4,np.nan)
        bammer5 = np.where(mask5 != 1.0, bammer5,np.nan)

        bammer1 = np.where(mask1,bammer1, np.nan)
        bammer2 = np.where(mask2,bammer2, np.nan)
        bammer3 = np.where(mask3,bammer3, np.nan)
        bammer4 = np.where(mask4,bammer4, np.nan)
        bammer5 = np.where(mask5,bammer5, np.nan)

        prop1 = np.where(mask1 != 1.0, prop1,np.nan)
        prop2 = np.where(mask2 != 1.0, prop2,np.nan)
        prop3 = np.where(mask3 != 1.0, prop3,np.nan)
        prop4 = np.where(mask4 != 1.0, prop4,np.nan)
        prop5 = np.where(mask5 != 1.0, prop5,np.nan)

        prop1 = np.where(mask1,prop1, np.nan)
        prop2 = np.where(mask2,prop2, np.nan)
        prop3 = np.where(mask3,prop3, np.nan)
        prop4 = np.where(mask4,prop4, np.nan)
        prop5 = np.where(mask5,prop5, np.nan)
    else:
        bammer1 = np.where(mask1 == 1.0, bammer1,np.nan) # csf
        bammer2 = np.where(mask2 == 1.0, bammer2,np.nan)
        bammer3 = np.where(mask3 == 1.0, bammer3,np.nan)
        bammer4 = np.where(mask4 == 1.0, bammer4,np.nan)
        bammer5 = np.where(mask5 == 1.0, bammer5,np.nan)

        bammer1 = np.where(mask1,bammer1, np.nan)
        bammer2 = np.where(mask2,bammer2, np.nan)
        bammer3 = np.where(mask3,bammer3, np.nan)
        bammer4 = np.where(mask4,bammer4, np.nan)
        bammer5 = np.where(mask5,bammer5, np.nan)

        prop1 = np.where(mask1 == 1.0, prop1,np.nan)
        prop2 = np.where(mask2 == 1.0, prop2,np.nan)
        prop3 = np.where(mask3 == 1.0, prop3,np.nan)
        prop4 = np.where(mask4 == 1.0, prop4,np.nan)
        prop5 = np.where(mask5 == 1.0, prop5,np.nan)

        prop1 = np.where(mask1,prop1, np.nan)
        prop2 = np.where(mask2,prop2, np.nan)
        prop3 = np.where(mask3,prop3, np.nan)
        prop4 = np.where(mask4,prop4, np.nan)
        prop5 = np.where(mask5,prop5, np.nan)
    
    return bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5
