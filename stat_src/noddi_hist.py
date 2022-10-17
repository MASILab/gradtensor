import numpy as np
from numpy import linspace
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import nibabel as nib
from noddi_utils import get_data

def plot_den(axes, bammer, prop, x, y):
    # plot histogram chart for var1
    prop = prop.flatten()
    bammer = bammer.flatten()
    bammer = bammer[~np.isnan(bammer)]
    prop = prop[~np.isnan(prop)]
    #sns.histplot(x=prop, bins=20, edgecolor='black', ax=axes[x,y])
    sns.histplot(x=prop, stat="density", bins=20, edgecolor='black', ax=axes[x,y])
    

    # plot histogram chart for var2
    n_bins = 20
    # get positions and heights of bars
    heights, bins = np.histogram(bammer, density=True, bins=n_bins) 
    # heights, bins = np.histogram(bammer, density=True, bins=n_bins) 

    # multiply by -1 to reverse it
    heights *= -1
    bin_width = np.diff(bins)[0]
    bin_pos =( bins[:-1] + bin_width / 2) * 1

    # plot
    axes[x,y].bar(bin_pos, heights, width=bin_width, edgecolor='black')


plt.rcParams.update({'font.size':15})
fig, axes = plt.subplots(5,3, figsize=(13,18))
s = ['ficvf.nii', 'fiso.nii', 'odi.nii']
for i in range(3):
    bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5 = get_data(s[i])

    plot_den(axes, bammer1, prop1,  0, i)
    plot_den(axes, bammer2, prop2,  1, i)
    plot_den(axes, bammer3, prop3,  2, i)
    plot_den(axes, bammer4, prop4,  3, i)
    plot_den(axes, bammer5, prop5,  4, i)

plt.subplots_adjust(wspace=0.5,hspace=0.23)
plt.show()
