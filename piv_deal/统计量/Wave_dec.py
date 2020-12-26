# 使用连续小波变换，将单点的脉动速度的时间序列进行尺度分解，并分为两个尺度
import numpy as np
import matplotlib.pyplot as plt
import pycwt as pt
import os
import time

nu = 1
NN = 8
N = NN * nu
m = 399
n = 65
fs = 200
dt = 1 / fs

fmat = '{}d'.format(m*n)
s = (N, m*n)
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
for i in range(0, N):
    path_fu = path_fuv + "fu_{:05}.bin".format(i+1)
    fu[i] = np.fromfile(path_fu, dtype=fmat)
end = time.perf_counter()
print(end)
print(end-start)
wave, scales, freqs, coi, fft, fftfreqs = pt.cwt(fu[1], dt, wavelet='morlet')
iwave = pt.icwt(wave, scales, dt, wavelet= 'morlet')
