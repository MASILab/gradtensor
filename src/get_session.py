import numpy as np
import os
import glob

p = []
for path in glob.glob('/nfs/masi/kanakap/projects/LR_tract/MASiVar_kids/*/*'):
    p.append(path)

d = {}
for i in range(len(p)):
    a = [p[i]+'/tracto_op_1_est',p[i]+'/tracto_op_2_est']
    d[i] = a


