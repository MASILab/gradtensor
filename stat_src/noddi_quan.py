import sys
import math
from scipy import stats
import numpy as np
import pandas as pd
import nibabel as nib
import seaborn as sns
from noddi_utils import get_data
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


def med_scalar(whole_brain):
    m = np.nanmedian(whole_brain)
    return m


def df_med(bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5):
    a = med_scalar(bammer1)
    b = med_scalar(bammer2)
    c = med_scalar(bammer3)
    d = med_scalar(bammer4)
    e = med_scalar(bammer5)

    f = med_scalar(prop1)
    g = med_scalar(prop2)
    h = med_scalar(prop3)
    i = med_scalar(prop4)
    j = med_scalar(prop5)

    bamm = [a, b, c, d, e]
    prop = [f, g, h, i, j]
    df = pd.DataFrame(bamm)
    df = df.assign(x = 'Empirical Bammer')
    dfp = pd.DataFrame(prop)
    dfp = dfp.assign(x = 'Proposed')

    cat = pd.concat([df,dfp])
    return cat, bamm, prop

def plot_sub(axes, cat, bamm, prop, n, title):#, ymin, ymax):

    palette= {'Empirical Bammer': 'teal','Proposed': 'rebeccapurple'}
    sns.swarmplot( ax = axes[n], x="x", y=0, hue="x", data=cat, palette=palette)
    idx0 = 0
    idx1 = 1
    locs1 = axes[n].get_children()[idx0].get_offsets()
    locs2 = axes[n].get_children()[idx1].get_offsets()
    sort_idxs1 = np.array([0, 1, 2, 3, 4])
    sort_idxs2 = np.array([0, 1, 2, 3, 4])#np.argsort(prop)
    locs2_sorted = locs2[sort_idxs2.argsort()][sort_idxs1]
    for i in range(locs1.shape[0]):
        x = [locs1[i, 0], locs2_sorted[i, 0]]
        y = [locs1[i, 1], locs2_sorted[i, 1]]
        axes[n].plot(x, y, '--', color="black", alpha=0.2)

    palette= {'Empirical Bammer': 'mediumturquoise','Proposed': 'mediumpurple'}

    bb = sns.boxplot( ax = axes[n], x="x", y=0, hue="x", data=cat,palette=palette, width = 0.4,dodge=False)

    # for i,artist in enumerate(bb.artists):
    # Set the linecolor on the artist to the facecolor, and set the facecolor to None
      #  col = artist.get_facecolor()
      #  artist.set_edgecolor(col)
      #  artist.set_facecolor('None')
      #  for j in range(i*6,i*6+6):
      #      line = bb.lines[j]
      #      line.set_color(col)
      #      line.set_mfc(col)
      #      line.set_mec(col)

    # axes[n].set(ylim=(ymin, ymax))
    axes[n].set_title(' ', fontsize=15)
    axes[n].set_ylabel('Median ' + title,fontsize = 15)
    axes[n].set_xlabel('')
    axes[n].set_xticklabels([])
    axes[n].tick_params(labelsize=15)
    axes[n].grid(True,linewidth=0.5)
    axes[n].legend([],[],frameon=False)


bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5 = get_data('ficvf.nii')
cat_iVF, bamm_iVF, prop_iVF = df_med(bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5)
bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5 = get_data('fiso.nii')
cat_cVF, bamm_cVF, prop_cVF  = df_med(bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5)
bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5 = get_data('odi.nii')
cat_ODI, bamm_ODI, prop_ODI  = df_med(bammer1,bammer2,bammer3,bammer4,bammer5,prop1,prop2,prop3,prop4,prop5)

fig, axes = plt.subplots(1,3, figsize=(15,5))
plot_sub(axes, cat_iVF, bamm_iVF, prop_iVF, 0, 'iVF')#,  0.078, 0.108)
plot_sub(axes, cat_cVF, bamm_cVF, prop_cVF, 1, 'cVF')#, 0.012, 0.025)
plot_sub(axes, cat_ODI, bamm_ODI, prop_ODI, 2, 'ODI')#, 0.073, 0.098)

axes[2].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
plt.subplots_adjust(wspace=0.35,hspace=0.1)
plt.show()

pivf = stats.wilcoxon(bamm_iVF, prop_iVF)
pcvf = stats.wilcoxon(bamm_cVF, prop_cVF)
podi = stats.wilcoxon(bamm_ODI, prop_ODI)
print(pivf,pcvf,podi)

