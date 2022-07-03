import re
import os
import glob 
import subprocess
p = []
for path in glob.glob('/nfs/masi/MASIver/data/production/MASiVar_PRODUCTION_v2.0.0/derivatives/prequal-v1.0.0/sub-cIVs*/*/dwi'):
    p.append(path)
p.sort()

count = 0
for i in p:
    print(i)
    run1 = {'dwi':[] ,'bvec': [], 'bval': []}
    run2 = {'dwi':[] ,'bvec': [], 'bval': []}
    for f in os.listdir(i):
        pattern = 'APP'
        if re.search(pattern,f):
            prun1 = 'run-1'
            if re.search(prun1,f):
                if f.endswith('.gz'):
                    l = i + "/" + f
                    run1['dwi'].append(l)
                if f.endswith('.bvec'):
                    l = i + "/" + f
                    run1['bvec'].append(l)
                if f.endswith('.bval'):
                    l = i + "/" + f
                    run1['bval'].append(l)
            prun2 = 'run-2'
            if re.search(prun2,f):
                if f.endswith('.gz'):
                    l = i + "/" + f
                    run2['dwi'].append(l)
                if f.endswith('.bvec'):
                    l = i + "/" + f
                    run2['bvec'].append(l)
                if f.endswith('.bval'):
                    l = i + "/" + f
                    run2['bval'].append(l)
    run1['dwi'].sort()
    run1['bvec'].sort()
    run1['bval'].sort()
    a = run1['dwi'][0]
    subj = a.split('/')[9]
    sess = a.split('/')[10]
    outdir = '/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/' + subj + "/"+ sess +'/prequal_dwi_cat/'
    out_dwi1 = outdir + subj + "_" + sess  + '_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.nii.gz'
    out_bval1 = outdir + subj + "_" + sess + '_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.bval'
    out_bvec1 = outdir + subj + "_" + sess + '_acq-b1000b2000n96r21x21x22peAPP_run-1_dwi.bvec'

    command1 = ["/home/local/VANDERBILT/kanakap/miniconda3/envs/dipyconda/bin/scil_concatenate_dwi.py", out_dwi1 , out_bval1 , out_bvec1 , "--in_dwis", run1['dwi'][0] , run1['dwi'][1] , "--in_bvals" , run1['bval'][0] , run1['bval'][1] , "--in_bvecs" , run1['bvec'][0] , run1['bvec'][1]]
    print(command1)
    proc = subprocess.Popen(command1)
    proc.wait()
    
    run2['dwi'].sort()
    run2['bvec'].sort()
    run2['bval'].sort()
    out_dwi2 = outdir + subj + "_" + sess + '_acq-b1000b2000n96r21x21x22peAPP_run-2_dwi.nii.gz'
    out_bval2 = outdir + subj + "_" + sess + '_acq-b1000b2000n96r21x21x22peAPP_run-2_dwi.bval'
    out_bvec2 = outdir + subj + "_" + sess + '_acq-b1000b2000n96r21x21x22peAPP_run-2_dwi.bvec'
    
    command2 = ["/home/local/VANDERBILT/kanakap/miniconda3/envs/dipyconda/bin/scil_concatenate_dwi.py",out_dwi2 , out_bval2 , out_bvec2 , "--in_dwis",run2['dwi'][0] , run2['dwi'][1] , "--in_bvals" , run2['bval'][0] ,  run2['bval'][1] ,"--in_bvecs",run2['bvec'][0] , run2['bvec'][1]]
    print(command2)
    proc2 = subprocess.Popen(command2)
    proc2.wait()

    count+=1
print(count)


