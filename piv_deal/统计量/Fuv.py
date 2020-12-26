# 此程序计算瞬时脉动速度文件
import numpy as np
import os
import time

nu = 2
NN = 8215
scale = 3.75
N = NN * nu
m = 399
n = 65
fmat = '{}d'.format(m*n)
s = m*n
fu = np.zeros(s)
fv = np.zeros(s)
path_data = "L:\\M_C\\1p\\"
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_Sta = path_result + 'Statistic\\'
path_mean = path_Sta + 'Velocity_mean.dat'
path_list = os.listdir(path_data1)
path_fuv = path_result + 'Fuv\\'
# 读取平均速度场
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
mean_x = mean_velocity[:, 0]
meany = mean_velocity[:, 1]
u = mean_velocity[:, 2]
v = mean_velocity[:, 3]
# 读取瞬时速度场
start = time.perf_counter()
for i in range(0, N):
    path_velocity = path_data1 + path_list[i]
    V = np.loadtxt(path_velocity, dtype='float', skiprows=3, comments='DATASETAUXDATA')
    fu = V[:, 2] * scale - u
    fv = V[:, 3] * scale - v
    path_fu = path_fuv + "fu_{:05}.bin".format(i+1)
    path_fv = path_fuv + "fv_{:05}.bin".format(i+1)
    fu.tofile(path_fu)
    fv.tofile(path_fv)
end = time.perf_counter()
