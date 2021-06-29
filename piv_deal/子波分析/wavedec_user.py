# 使用小波变换得到小波系数（空间）
import numpy as np
import time
import os
from Read_Parameter import read_p

if __name__ == "__main__":
    read_p
###################################################################
# read parameters from file
path_data = "L:\\M_C\\1p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8  # number of parameters
nu = 1  # number of cases
bu = read_p(para_metre, p_n)
mul = float(bu[0])  # coefficient of kinematic viscosity
dyn = int(bu[1])  # Deviance from zero to real wall coordinate
ps = int(bu[2])  # first point in log law region
pe = int(bu[3])  # last point in log law region
NN = int(bu[4])  # number of samples
m = int(bu[5])  # number of stream-wise points
n = int(bu[6])  # number of wall-normal points
scale = float(bu[7])  # amplification coefficient
##################################################################
# basic parameter
N = NN * nu  # sum of all frames
Mx = 8  # 分解最大尺度
fs = 200  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
h = [0.125, 0.375, 0.375, 0.125]
g = [-0.00610, -0.0869, -0.5789, -0.5789, 0.0869, 0.00610]
format_fuv = '{}d'.format(m * n)
s = (N, m * n)
# fv = np.zeros(s)
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_Sta = path_result + 'Statistic\\'
path_list = os.listdir(path_data1)
path_fuv = path_result + 'Fuv\\'
path_mean = path_Sta + 'Velocity_mean.dat'
path_stat = path_Sta + 'information.dat'
PSD = path_Sta + 'PSD//'
wave_e = path_result + 'coef\\'
wave_f = path_result + 'wavet\\'
yu = read_p(path_stat, 5)
uto = yu[3]
bc_th = yu[0]

start = time.perf_counter()
print(start)


def stay(dx, ik):
    if ik < 0:
        ln = -ik - 1
    else:
        if ik > dx - 1:
            ln = 2 * dx - 1 - ik
        else:
            ln = ik
    return ln


# 读取平均速度场
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meany = mean_velocity[:, 1].reshape(m, n, order='F')[1, :]
x = mean_velocity[:, 0].reshape(m, n, order='F')
y = mean_velocity[:, 1].reshape(m, n, order='F')
##################################################################
# read velocity fields
eu = np.zeros((n + 1, m + 1))
ev = np.zeros((n + 1, m + 1))
for ii in range(0, N):
    path_fu = path_fuv + "fu_{:05}.bin".format(ii + 1)
    path_fv = path_fuv + "fv_{:05}.bin".format(ii + 1)
    fu = np.fromfile(path_fu, dtype=format_fuv)
    fv = np.fromfile(path_fv, dtype=format_fuv)
    fuu = fu.reshape(m, n, order='F')
    fvv = fv.reshape(m, n, order='F')
    e1 = np.zeros((n, m + 1))
    e2 = np.zeros((n, m + 1))
    c = 1
    for k in range(Mx):
        b1 = np.zeros((m, n))
        b2 = np.zeros((m, n))
        c1 = np.zeros((m, n))
        c2 = np.zeros((m, n))
        for iy in range(n):
            for i in range(m):
                for j in range(-1, 3):
                    b1[i, iy] = b1[i, iy] + h[j + 1] * fuu[stay(m, i - c * j), iy]
                    b2[i, iy] = b2[i, iy] + h[j + 1] * fvv[stay(m, i - c * j), iy]
                for j in range(-2, 4):
                    c1[i, iy] = c1[i, iy] + g[j + 2] * fuu[stay(m, i - c * j), iy]
                    c2[i, iy] = c2[i, iy] + g[j + 2] * fvv[stay(m, i - c * j), iy]
                e1[iy, k] = e1[iy, k] + c1[i, iy] * c1[i, iy]
                e2[iy, k] = e2[iy, k] + c2[i, iy] * c2[i, iy]
            eu[iy, k] = eu[iy, k] + e1[iy, k] / m
            ev[iy, k] = ev[iy, k] + e2[iy, k] / m
        fuu = b1
        fvv = b2
        c = c * 2
        cc1 = c1.reshape(m*n, 1, order='F')
        cc2 = c2.reshape(m*n, 1, order='F')
        path_wave_eu = wave_e + "Nu{:04}_sacle{:02}.bin".format(ii, k)
        path_wave_ev = wave_e + "Nv{:04}_sacle{:02}.bin".format(ii, k)
        cc1.tofile(path_wave_eu)
        cc2.tofile(path_wave_ev)

euu = np.sum(eu, 1)
evv = np.sum(ev, 1)
path_wave_e = wave_f + "energy.dat"
f = open(path_wave_e, 'w')
f.write('VARIABLES= "ix", "iy", "cu", "cv" \n')
f.write('Zone T="%s" i=%d j= %d\n' % ('1pumps', Mx, n))
for iy in range(n):
    for k in range(Mx):
        f.write("{}\t{}\t{}\t{}\n".format(k, iy, eu[iy, k] / euu[iy], ev[iy, k] / evv[iy]))
f.close()
