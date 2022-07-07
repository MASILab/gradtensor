import numpy as np
import os
import glob
import json 

p = []
for path in glob.glob('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/*/*'):
    p.append(path)

d = {}
for i in range(len(p)):
    a = [p[i]+'/tracto_op_1_est/stats_json/length_stats.json',p[i]+'/tracto_op_2_est/stats_json/length_stats.json']
    d[i] = a

all_len = = {}

for i in d:
    if os.path.exists(d[i][0]):
        f = open(d[i][0])
        data = json.load(f)

    if os.path.exists(d[i][1]):
        f = open(d[i][1])
        data1 = json.load(f)
    
    k = list(data.keys())[0]
    k1 = list(data1.keys())[0]
    if k == k1:
        v = list(data[k].values())
        v1 = list(data1[k1].values())
        su = 0
        if len(v) == len(v1):
            for j in range(len(v)):
                diff = v[j]['mean_length'] - v1[j]['mean_length']
                su += diff
            mean_diff = su / len(v)

    all_len[i] = [mean_diff]

print(all_len)

#~/MATLAB_2017b_install/bin/matlab
        



