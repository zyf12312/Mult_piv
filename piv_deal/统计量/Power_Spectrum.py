# calculate power spectral density for piv data
import numpy as np
from scipy.fftpack import fft
import math
import os
import time
from Read_Parameter import read_p
if __name__ == "__main__":
    read_p
###################################################################
# read parameters from file
path_data = "L:\\M_C\\4p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8                    # number of parameters
nu = 1                     # number of cases
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
fs = 800  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
##################################################################
# Set path
format_fuv = '{}d'.format(m * n)
s = (N, m * n)
fu = np.zeros(s)
# fv = np.zeros(s)
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_Sta = path_result + 'Statistic\\'
path_list = os.listdir(path_data1)
path_fuv = path_result + 'Fuv\\'
path_mean = path_Sta + 'Velocity_mean.dat'
path_stat = path_Sta + 'information.dat'
PSD = path_Sta + 'PSD//'
yu = read_p(path_stat, 5)
uto = yu[3]
bc_th = yu[0]
# 读取平均速度场
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meany = mean_velocity[:, 1].reshape(m, n, order='F')[1, :]
u = mean_velocity[:, 2].reshape(m, n, order='F')
um = np.mean(u, 0)
# 读取壁面摩擦速度

start = time.perf_counter()
print(start)
##################################################################
# read files
for i in range(0, N):
    path_fu = path_fuv + "fu_{:05}.bin".format(i + 1)
    fu[i] = np.fromfile(path_fu, dtype=format_fuv)
end = time.perf_counter()
print(end)
print(end - start)
f1 = np.arange(N)
f = f1[range(int(N/2))] / N * fs
for i in range(n):
    abs_fft = np.zeros(int(N/2))
    for j in range(m):
        tf = fft(fu[:, i*m+j])
        abs_ft = np.abs(tf)**2 / N / N
        abs_fft = abs_fft + abs_ft[range(int(N/2))]/m
    k = 2 * math.pi * f / um[i] * bc_th
    abs_fft = abs_fft / uto / uto
    aa = abs_fft * k
    ft = open(PSD + "PSD{:02}.dat".format(i + 1), 'w')
    ft.write('VARIABLES= "k<sub>u</sub><greek>d</greek>", "<greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek><sup>2</sup>", "k<sub>u</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek><sup>2</sup>"\n')
    ft.write('Zone T="%s" i=%d\n' % ('1pumps', int(N/2)))
    for j in range(0, int(N/2)):
        ft.write("{}\t{}\t{}\n".format(k[j], abs_fft[j], aa[j]))
    ft.close()
