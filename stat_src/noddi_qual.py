import sys
import math
import numpy as np
import nibabel as nib
import seaborn as sns
from noddi_utils import get_data
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_slice(image,i):#, mask):
    slice_idx = 60
    #image = np.where(mask, image,0)
    slice = image[slice_idx, :,:]
    slice = np.flip(np.rot90(slice,3))
    slice = np.nan_to_num(slice)
    plt.subplot(3,5,i)
    plt.axis('off')
    cmap = plt.get_cmap('gray')
    m = 0
    M = 1
    im = plt.imshow(slice, vmin=m, vmax=M, cmap=cmap)

fig, axes = plt.subplots(3,5, figsize=(16,5))
plt.rcParams.update({'font.size':15})
bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5 = get_data(sys.argv[1])

d = [bammer1,bammer2,bammer3,bammer4,bammer5,
    prop1,prop2,prop3,prop4,prop5,
    np.abs(bammer1 - prop1), np.abs(bammer2 - prop2), np.abs(bammer3 - prop3), np.abs(bammer4 - prop4),  np.abs(bammer5 - prop5)]

for i in range(15):
    plot_slice(d[i], i+1)

plt.subplots_adjust(wspace=0.01,hspace=0.001) 
plt.show()
