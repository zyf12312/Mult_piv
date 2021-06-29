# 使用连续小波变换，将单点的脉动速度的时间序列进行尺度分解，并分为两个尺度
import numpy as np
import matplotlib.pyplot as plt
import pycwt as pt
import os
import time
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
##################################################################
# basic parameter

N = NN * nu  # sum of all frames
fs = 200  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
##################################################################
# Set path
format_fuv = '{}d'.format(m * n)
s = (N, m * n)
fu = np.zeros(s)
# fv = np.zeros(s)
path_data = "L:\\M_C\\1p\\"
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_Sta = path_result + 'Statistic\\'
path_mean = path_Sta + 'Velocity_mean.dat'
path_list = os.listdir(path_data1)
path_fuv = path_result + 'Fuv\\'
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
##################################################################
# wavelet transform
wave, scales, freqs, coi, fft, fftfreqs = pt.cwt(fu[:, 0], dt, 0.1, dt * 2, wavelet='morlet')
iwave = pt.icwt(wave, scales, dt, 0.1, wavelet='morlet')
plt.plot(fu[:, 0])
plt.plot(iwave.real)
plt.show
