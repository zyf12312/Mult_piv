# 此程序求取雷诺应力及脉动速度场
import numpy as np
import os
from Read_Parameter import read_p
if __name__ == "__main__":
    read_p
path_data = "L:\\M_C\\1p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8
nu = 1
bu = read_p(para_metre, p_n)
mul = float(bu[0])         # coefficient of kinematic viscosity
dyn = int(bu[1])           # Deviance from zero to real wall coordinate
ps = int(bu[2])            # first point in log law region
pe = int(bu[3])            # last point in log law region
NN = int(bu[4])            # number of samples
m = int(bu[5])             # number of stream-wise points
n = int(bu[6])             # number of wall-normal points
scale = float(bu[7])       # amplification coefficient
s = m*n
N = NN * nu
fuu = np.zeros(s)
fuv = np.zeros(s)
fvv = np.zeros(s)
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\Statistic\\"
path_mean = path_result + 'Velocity_mean.dat'
path_list = os.listdir(path_data1)
path_tau = path_result + 'Mu_tau.dat'
path_Reynold_Stress = path_result + 'Reynold_Stress.dat'
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
mean_x = mean_velocity[:, 0]
meany = mean_velocity[:, 1]
u = mean_velocity[:, 2]
v = mean_velocity[:, 3]
um = np.mean(u, 0)
for i in range(0, N):
    path_velocity = path_data1 + path_list[i]
    V = np.loadtxt(path_velocity, dtype='float', skiprows=3, comments='DATASETAUXDATA')
    fu = V[:, 2] * scale - u
    fv = V[:, 3] * scale - v
    fuu = fuu + fu**2
    fvv = fvv + fv**2
    fuv = fuv + fu*fv
fuu = fuu / N
fuv = fuv / N
fvv = fvv / N
f = open(path_Reynold_Stress, 'w')
f.write('VARIABLES= "x(mm)", "y(mm)", "uu", "vv", "uv" \n')
f.write('Zone T="%s" i=%d j= %d\n' % ('1pumps', m, n))
for i in range(0, m*n):
    f.write("{}\t{}\t{}\t{}\t{}\n".format(mean_x[i], meany[i], fuu[i], fvv[i], fuv[i]))
f.close()
