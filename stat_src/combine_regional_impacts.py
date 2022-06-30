import pandas as pd
import glob 

ip_dir = "/nfs/masi/kanakap/projects/LR/masivar_output/SNRinf_d32_1"

def read_rename_csv(filename, matter, metric):
    csvfile = pd.read_csv(ip_dir + '/' + filename)
    csvfile = csvfile.rename(columns={'Unnamed: 0': matter, '0': metric})
    return csvfile

csv1 = read_rename_csv('fa_gm.csv','GM','FA')
csv2 = read_rename_csv('md_gm.csv','GM','MD')
csv3 = read_rename_csv('rd_gm.csv','GM','RD')
csv4 = read_rename_csv('ad_gm.csv','GM','AD')
csv5 = read_rename_csv('v1_gm.csv','GM','V1')


csv6 = read_rename_csv('fa_wm.csv','WM','FA')
csv7 = read_rename_csv('md_wm.csv','WM','MD')
csv8 = read_rename_csv('rd_wm.csv','WM','RD')
csv9 = read_rename_csv('ad_wm.csv','WM','AD')
csv10 = read_rename_csv('v1_wm.csv','WM','V1')


gm_merged = pd.concat([csv1,csv2,csv3,csv4,csv5],axis=1)
gm_merged = gm_merged.loc[:,~gm_merged.columns.duplicated()].copy()
wm_merged = pd.concat([csv6,csv7,csv8,csv9,csv10],axis=1)
wm_merged = wm_merged.loc[:,~wm_merged.columns.duplicated()].copy()

gm_merged.to_csv('new_lr_impact_gm.csv')
wm_merged.to_csv('new_lr_impact_wm.csv')




