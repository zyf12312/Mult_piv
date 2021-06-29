# Decomposition velocities fields by Proper Orthogonal Decomposition
import numpy as np
import time
from mrpod import pod_modes
from Read_Parameter import read_p
if __name__ == "__main__":
    read_p
###################################################################
# read parameters from file
path_data = "L:\\M_C\\4p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8                    # number of parameters
nu = 2                     # number of cases
bu = read_p(para_metre, p_n)
mul = float(bu[0])         # coefficient of kinematic viscosity
dyn = int(bu[1])           # Deviance from zero to real wall coordinate
ps = int(bu[2])            # first point in log law region
pe = int(bu[3])            # last point in log law region
NN = int(bu[4])            # number of samples
m = int(bu[5])             # number of stream-wise points
n = int(bu[6])             # number of wall-normal points
scale = float(bu[7])       # amplification coefficient
##################################################################
# basic parameter
N = NN * nu  # sum of all frames
N = 100
fs = 200  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
nm = 50
start = time.perf_counter()
##################################################################
# Set path
format_fuv = '{}d'.format(m * n)
s = (N, 2 * m * n)
fu = np.zeros(s)
path_result = path_data + "Result\\"
path_fuv = path_result + 'Fuv\\'
path_pod = path_result + 'POD\\'

path_mean = path_result + 'Statistic\\Velocity_mean.dat'
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meanx = mean_velocity[:, 0]
meany = mean_velocity[:, 1]
for i in range(0, N):
    path_fu = path_fuv + "fu_{:05}.bin".format(i + 1)
    path_fv = path_fuv + "fv_{:05}.bin".format(i + 1)
    fu[i, 0:m*n] = np.fromfile(path_fu, dtype=format_fuv)
    fu[i, m*n:2*m*n] = np.fromfile(path_fv, dtype=format_fuv)
print("开始计算")
pod_result = pod_modes(fu, num_of_modes=nm, normalize_mode=True)
proj_co = pod_result['proj_coeffs']
mode = pod_result['modes']
eigv = pod_result['eigvals']
eiga = eigv/eigv.sum()*100
print("计算结束")
# for i in range(nm):
#     path_mode = path_pod + "2mod{:02}.dat".format(i)
#     fip = open(path_mode, 'w')
#     fip.write('VARIABLES= "x(mm)","y(mm)","u(m/s)","v(m/s)"\n')
#     fip.write('Zone T="%d" i=%d j=%d\n' % (i, m, n))
#     for j in range(m*n):
#         fip.write("{}\t{}\t{}\t{}\n".format(meanx[j], meany[j], mode[i, j], mode[i, m*n+j]))
#     fip.close()
# path_coer = path_pod + "6convergence.dat"
# fid = open(path_coer, 'w')
# fid.write('VARIABLES="i_level","energy(%)"\n')
# fid.write('Zone T="1" i=%d\n' % N)
# for i in range(N):
#     fid.write("{}\t{}\n".format(i, eiga[i]))
# fid.close()
# end_time = time.perf_counter()
# print(end_time-start)
