import matplotlib.pyplot as plt
import seaborn as sn
import scipy.io

mat = scipy.io.loadmat('/home/local/VANDERBILT/kanakap/gradtensor_var/matlab.mat')

theta = mat.get('theta')
theta_all = mat.get('theta_all')
FA = mat.get('FA')
FA_all = mat.get('FA_all')
Ldet = mat.get('Ldet')
Ldet_all = mat.get('Ldet_all')


FA_sim_corpt_x = mat.get('FA_sim_corpt_x')
FA_diff = mat.get('baba')


theta_all[theta_all < 18] = 0
theta_all[ (theta_all >= 18*1) & (theta_all < 18*2) ] = 1
theta_all[ (theta_all >= 18*2) & (theta_all < 18*3) ] = 2
theta_all[ (theta_all >= 18*3) & (theta_all < 18*4) ] = 3
theta_all[ (theta_all >= 18*4) & (theta_all < 18*5) ] = 4
theta_all[ (theta_all >= 18*5) & (theta_all < 18*6) ] = 5
theta_all[ (theta_all >= 18*6) & (theta_all < 18*7) ] = 6
theta_all[ (theta_all >= 18*7) & (theta_all < 18*8) ] = 7
theta_all[ (theta_all >= 18*8) & (theta_all < 18*9) ] = 8
theta_all[ (theta_all >= 18*9) & (theta_all < 18*10) ] = 9

colors = {'0':'red', '2':'green', '3':'blue', '4':'yellow'}
plt.scatter(FA_sim_corpt_x, FA_diff, c = theta_all)
plt.legend()

plt.show()

exit()
for idx, val in enumerate(theta_all):
    if val < 18*1:
        theta_all[idx,:] = 0
    if 18*1 <= val < 18*2:
        theta_all[idx,:] = 1
    if 18*2 <= val < 18*3:
        theta_all[idx,:] = 2
    if 18*3 <= val < 18*4:
        theta_all[idx,:] = 3
    if 18*4 <= val < 18*5:
        theta_all[idx,:] = 4
    if 18*5 <= val < 18*6:
        theta_all[idx,:] = 5
    if 18*6 <= val < 18*7:
        theta_all[idx,:] = 6
    if 18*7 <= val < 18*8:
        theta_all[idx,:] = 7
    if 18*8 <= val < 18*9:
        theta_all[idx,:] = 8
    if 18*9 <= val < 18*10:
        theta_all[idx,:] = 9

plt.scatter(FA_sim_corpt_x, FA_diff, c = theta_all)

plt.show()

