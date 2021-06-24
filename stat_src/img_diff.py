import nibabel as nib
import numpy as np

img = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/OUTPUTS_future_fieldmap/p_3tb_posA_mask_md.nii')
mask = nib.load('/nfs/masi/kanakap/gradtensor_data/PVP_phantom/scans/3tb/posA/mask.nii')

img_arr = np.array(img.dataobj)
mask_arr = np.array(mask.dataobj)
unique(mask_arr)
#print(len(img_arr))
#print(len(mask_arr))
#print(mask_arr)

#data = img_arr[mask_arr] #get only values in mask
#variance = np.var(data) #variance
#std = np.std(data) #standard deviation
#mean = np.mean(data) #mean (edited) 
