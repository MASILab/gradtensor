import nibabel as nib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
import pandas as pd
import math

true_fa =  nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_d32_/true__fa.nii').get_fdata()
noise_fainf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
corpt_fainf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1_corr/uncorrected_fa.nii').get_fdata()
bx_corr_fainf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fainf = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNRinf_d32_1_corr/approx_corrected_fa.nii').get_fdata()

noise_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
corpt_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1_corr/uncorrected_fa.nii').get_fdata()
bx_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa100 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR100_d32_1_corr/approx_corrected_fa.nii').get_fdata()

noise_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
corpt_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/uncorrected_fa.nii').get_fdata()
bx_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa30 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR30_d32_1_corr/approx_corrected_fa.nii').get_fdata()

noise_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1_corr/uncorrected_Nest_fa.nii').get_fdata()
corpt_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1_corr/uncorrected_fa.nii').get_fdata()
bx_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1_corr/emp/emp_corrected_fa.nii').get_fdata()
sm_corr_fa10 = nib.load('/nfs/masi/kanakap/projects/LR/aggregate_study/OUTPUT_masivar_SNR10_d32_1_corr/approx_corrected_fa.nii').get_fdata()
mask = nib.load('/home/local/VANDERBILT/kanakap/MASIVAR_LR_input1/mask.nii').get_fdata()

def offset_violinplot_halves(ax, delta, width, inner, direction):
    # offset stuff
    if inner == 'sticks':
        lines = ax.get_lines()
        for line in lines:
            if direction == 'horizontal':
                data = line.get_ydata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is top, move neg, direction backwards for horizontal
                    data -= delta
                else:
                    # type is bottom, move pos, direction backward for hori
                    data += delta
                line.set_ydata(data)
            elif direction == 'vertical':
                data = line.get_xdata()
                print(data)
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is left, move neg
                    data -= delta
                else:
                    # type is left, move pos
                    data += delta
                line.set_xdata(data)


    for ii, item in enumerate(ax.collections):
        # axis contains PolyCollections and PathCollections
        if isinstance(item, matplotlib.collections.PolyCollection):
            # get path
            path, = item.get_paths()
            vertices = path.vertices
            half_type = _wedge_dir(vertices, direction)
            # shift x-coordinates of path
            if half_type in ['top','bottom']:
               if inner in ["sticks", None]:
                    if half_type == 'top': # -> up
                        vertices[:,1] -= delta
                    elif half_type == 'bottom': # -> down
                        vertices[:,1] += delta
            elif half_type in ['left', 'right']:
                if inner in ["sticks", None]:
                    if half_type == 'left': # -> left
                        vertices[:,0] -= delta
                    elif half_type == 'right': # -> down
                        vertices[:,0] += delta

def _wedge_dir(vertices, direction):
    if direction == 'horizontal':
        result = (direction, len(set(vertices[1:5,1])) == 1)
    elif direction == 'vertical':
        result = (direction, len(set(vertices[-3:-1,0])) == 1)
    outcome_key = {('horizontal', True): 'bottom',
                   ('horizontal', False): 'top',
                   ('vertical', True): 'left',
                   ('vertical', False): 'right'}
    # if the first couple x/y values after the start are the same, it
    #  is the input direction. If not, it is the opposite
    return outcome_key[result]

def plot_one_snr(corpt_fa,true_fa,label,x):
    err_fa = corpt_fa - true_fa
    pe_fa = (err_fa / true_fa) * 100
    ape_fa = np.abs(pe_fa)
    ape_fa[np.isnan(ape_fa)] = 0
    ape_fa = np.where(mask,ape_fa,math.nan)
    ape_fa = ape_fa.reshape(-1)
    ape_fa = [x for x in ape_fa if math.isnan(x) == False]
    df_ape_fa = pd.DataFrame(ape_fa).assign(label = label)
    #df_ape_fa = df_ape_fa.assign(x = x)
    df_ape_fa = df_ape_fa.rename(columns={0:x})
    return df_ape_fa

lr_noiseinf = plot_one_snr(corpt_fainf,true_fa,'Uncorrected','emp_inf')
corr_lrinf = plot_one_snr(bx_corr_fainf,noise_fainf,'Corrected','emp_inf')

lr_noise100 = plot_one_snr(corpt_fa100,true_fa,'Uncorrected','emp_100')
corr_lr100 = plot_one_snr(bx_corr_fa100,noise_fa100,'Corrected','emp_100')

lr_noise30 = plot_one_snr(corpt_fa30,true_fa,'Uncorrected','emp_30')
corr_lr30 = plot_one_snr(bx_corr_fa30,noise_fa30,'Corrected','emp_30')

lr_noise10 = plot_one_snr(corpt_fa10,true_fa,'Uncorrected','emp_10')
corr_lr10 = plot_one_snr(bx_corr_fa10,noise_fa10,'Corrected','emp_10')

alr_noiseinf = plot_one_snr(corpt_fainf,true_fa,'Uncorrected','app_inf')
acorr_lrinf = plot_one_snr(sm_corr_fainf,noise_fainf,'Corrected','app_inf')

alr_noise100 = plot_one_snr(corpt_fa100,true_fa,'Uncorrected','app_100')
acorr_lr100 = plot_one_snr(sm_corr_fa100,noise_fa100,'Corrected','app_100')

alr_noise30 = plot_one_snr(corpt_fa30,true_fa,'Uncorrected','app_30')
acorr_lr30 = plot_one_snr(sm_corr_fa30,noise_fa30,'Corrected','app_30')

alr_noise10 = plot_one_snr(corpt_fa10,true_fa,'Uncorrected','app_10')
acorr_lr10 = plot_one_snr(sm_corr_fa10,noise_fa10,'Corrected','app_10')


fig, ax = plt.subplots(1,1)
sns.set(font_scale = 2)
#pd_data = pd.concat([lr_noiseinf, corr_lrinf, lr_noise100, corr_lr100, lr_noise30, corr_lr30, lr_noise10, corr_lr10])
#pd_data = pd_data.rename(columns={0:'Abs percent error (%)'})
pd_data = pd.concat([lr_noiseinf, lr_noise100, lr_noise30, lr_noise10, corr_lrinf, corr_lr100, corr_lr30, corr_lr10, alr_noiseinf, alr_noise100, alr_noise30, alr_noise10, acorr_lrinf, acorr_lr100, acorr_lr30, acorr_lr10])
pd_data = pd_data.melt(id_vars=['label'], value_vars=['emp_inf','emp_100','emp_30','emp_10','app_inf','app_100','app_30','app_10'],var_name='H', value_name='value')
ax = sns.violinplot(data=pd_data, hue = 'label', x = 'H',y='value',gridsize=2500,split=True,dodge=False,inner=None,linewidth=2,scale="width")
inner=None
# offset stuff
delta = 0.02
for ii, item in enumerate(ax.collections):
    # axis contains PolyCollections and PathCollections
    if isinstance(item, matplotlib.collections.PolyCollection):
        # get path
        path, = item.get_paths()
        vertices = path.vertices

        # shift x-coordinates of path
        if not inner:
            if ii % 2: # -> to right
                vertices[:,0] += delta
            else: # -> to left
                vertices[:,0] -= delta
        else: # inner='box' adds another type of PollyCollection
            if ii % 3 == 0:
                vertices[:,0] -= delta
            elif ii % 3 == 1:
                vertices[:,0] += delta
            else: # ii % 3 = 2
                pass
#sns.despine(left=True)
#ax.set(xlabel='SNR', ylabel = 'Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )')
ax.set_xlabel('SNR',fontsize = 17)
ax.set_ylabel('Abs percent error in FA (%) abs ( ( a - b / b ) * 100 )',fontsize = 17)
ax.set_xticklabels(['Inf','100','30','10','Inf','100','30','10'])
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.legend(loc='upper left')
plt.title('Comparison of correction technique on FA in whole brain')
plt.ylim([-1,4])
plt.tight_layout()
plt.grid()
plt.show()
