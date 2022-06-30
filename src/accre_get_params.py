import sys
import re
import os
from glob import glob
import subprocess

def get_params():
    d = glob("/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/*/*/prequal_dwi_cat/*", recursive = True)
    d.sort()
    p = []
    for i in d:
        j = i.replace('.','/')
        la = j.split('/')
        if la[11] == 'nii':
            if 'APP' in la[10]:
                p.append(la)
    final = ()
    for i in p:
        print("************************************************************")
        #print(i)
        mbvalue = re.search('acq-b(.+?)n',i[10])
        bvalue = mbvalue.group(1)
        mrun = re.search('run-(.+?)_dwi',i[10])
        run = mrun.group(1)
        sub_sess = i[7]+'/'+i[8]
        sub = i[7]
        sess = i[8]
        dwi = i[10] + '.nii'
        bvec = i[10] + '.bvec'
        bval = i[10] + '.bval'
        a = (sub_sess, dwi, bval, bvec, run, sub, sess)
        name = 'accre_slurm_scripts/'+ sub + "_" + sess + "_" + bvalue + "_" + run +'.txt'
        with open(name, 'w') as f: 
            f.write("/nobackup/p_masi_brain_map/kanakap/LR/scripts/final_pipeline_LR_prequal.sh " + sub_sess +" "+ dwi +" "+ bval +" "+ bvec +" "+ run +" "+ sub +" " + sess)
            f.write("\n")
        #print("/scratch/kanakap/LR/log/"+ sub + "_" + sess + "_" + bvalue + "_" + run + ".out")
        #print("/nobackup/p_masi_brain_map/kanakap/LR/scripts/final_pipeline_LR_prequal.sh " + sub_sess +" "+ dwi +" "+ bval +" "+ bvec +" "+ run +" "+ sub +" " + sess)
        #proc = subprocess.Popen(["bash","/home/local/VANDERBILT/kanakap/gradtensor/src/final_pipeline_LR_prequal.sh", sub_sess, dwi, bval, bvec, bvalue, run, sub, sess]) 
        #proc.wait()
get_params()




