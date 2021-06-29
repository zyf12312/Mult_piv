# calculate power spectral density for piv data
import numpy as np
import os
import time
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
fs = 200  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
##################################################################
# Set path
format_fuv = '{}d'.format(m * n)
s = (N, m * n)
Mx = 8  # 分解最大尺度
c = 16
cy = 8
Nx = 2 * c
Ny = 2 * cy
fu = np.zeros(s)
fv = np.zeros(s)
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_Sta = path_result + 'Statistic\\'
path_list = os.listdir(path_data1)
path_fuv = path_result + 'Fuv\\'
path_mean = path_Sta + 'Velocity_mean.dat'
path_stat = path_Sta + 'information.dat'
wave_e = path_result + 'coef\\'
wave_f = path_result + 'wavet\\'
yu = read_p(path_stat, 5)
uto = yu[1]
bc_th = yu[0]
# 读取平均速度场
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meany = mean_velocity[:, 1].reshape(m, n, order='F')[1, :]
x = mean_velocity[:, 0].reshape(m, n, order='F')
xm = np.mean(x, 1)
dx = (xm[2] - xm[1]) / 1000
# 读取壁面摩擦速度

start = time.perf_counter()
print(start)
##################################################################
# read fluctuate velocities fields
for i in range(0, N):
    path_fu = path_fuv + "fu_{:05}.bin".format(i + 1)
    path_fv = path_fuv + "fv_{:05}.bin".format(i + 1)
    fu[i] = np.fromfile(path_fu, dtype=format_fuv)
    fv[i] = np.fromfile(path_fv, dtype=format_fuv)


for i in range(Mx):  # 尺度循环
    fu2p = np.zeros((n-Ny, Nx, Ny))
    fv2p = np.zeros((n-Ny, Nx, Ny))
    fuup = np.zeros((n-Ny, Nx, Ny))
    fvvp = np.zeros((n-Ny, Nx, Ny))
    fuvp = np.zeros((n-Ny, Nx, Ny))
    fu2n = np.zeros((n-Ny, Nx, Ny))
    fv2n = np.zeros((n-Ny, Nx, Ny))
    fuun = np.zeros((n-Ny, Nx, Ny))
    fvvn = np.zeros((n-Ny, Nx, Ny))
    fuvn = np.zeros((n-Ny, Nx, Ny))
    p1 = np.zeros((n-Ny, 1))
    n1 = np.zeros((n-Ny, 1))
    for j in range(N):  # 时间循环
        fu1 = fu[j].reshape(m, n, order='F')
        fv1 = fv[j].reshape(m, n, order='F')
        path_wave_eu = wave_e + "Nu{:04}_sacle{:02}.bin".format(j, i)
        wu = np.fromfile(path_wave_eu, dtype=format_fuv).reshape(m, n, order='F')
        for k in range(0, n - Ny):  # 法向循环
            for kk in range(c, m - c):  # 流向循环
                if (wu[kk, k] > 0) & (wu[kk, k] > wu[kk + 1, k]) & (wu[kk, k] > wu[kk - 1, k]):
                    fu2p[k, :, :] += fu1[kk - c:kk + c, k:k + Ny]
                    fv2p[k, :, :] += fv1[kk - c:kk + c, k:k + Ny]
                    fuup[k, :, :] += fu1[kk - c:kk + c, k:k + Ny] * fu1[kk - c:kk + c, k:k + Ny]
                    fvvp[k, :, :] += fv1[kk - c:kk + c, k:k + Ny] * fv1[kk - c:kk + c, k:k + Ny]
                    fuvp[k, :, :] += fu1[kk - c:kk + c, k:k + Ny] * fv1[kk - c:kk + c, k:k + Ny]
                    p1[k] += 1
                elif (wu[kk, k] < 0) & (wu[kk, k] < wu[kk + 1, k]) & (wu[kk, k] < wu[kk - 1, k]):
                    fu2n[k, :, :] += fu1[kk - c:kk + c, k:k + Ny]
                    fv2n[k, :, :] += fv1[kk - c:kk + c, k:k + Ny]
                    fuun[k, :, :] += fu1[kk - c:kk + c, k:k + Ny] * fu1[kk - c:kk + c, k:k + Ny]
                    fvvn[k, :, :] += fv1[kk - c:kk + c, k:k + Ny] * fv1[kk - c:kk + c, k:k + Ny]
                    fuvn[k, :, :] += fu1[kk - c:kk + c, k:k + Ny] * fv1[kk - c:kk + c, k:k + Ny]
                    n1[k] += 1
    for k in range(n-Ny):
        if p1[k] > 0:
            fu2p[k, :, :] /= p1[k] * uto
            fv2p[k, :, :] /= p1[k] * uto
            fuup[k, :, :] /= p1[k] * uto * uto
            fuvp[k, :, :] /= p1[k] * uto * uto
            fvvp[k, :, :] /= p1[k] * uto * uto
        if n1[k] > 0:
            fu2n[k, :, :] /= n1[k] * uto
            fv2n[k, :, :] /= n1[k] * uto
            fuun[k, :, :] /= n1[k] * uto * uto
            fuvn[k, :, :] /= n1[k] * uto * uto
            fvvn[k, :, :] /= n1[k] * uto * uto
        path_wave_e = wave_f + "Scale{:04}_layer{:02}.dat".format(i, k)
        f = open(path_wave_e, 'w')
        f.write('VARIABLES= "dx/D", "dy/D", "fu_p/u", "fv_p/u" "fuu_p/uu", "fuv_p/uu", "fvv_p/uu", "fu_n" "fv_n", "fuu_n", "fuv_n", "fvv_n" \n')
        f.write('Zone T="%s" i=%d j= %d\n' % ('1pumps', Nx, Ny))
        for iy in range(Ny):
            for kk in range(Nx):
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(dx*kk / bc_th, dx*iy / bc_th, fu2p[k, kk, iy], fv2p[k, kk, iy], fuup[k, kk, iy], fuvp[k, kk, iy], fvvp[k, kk, iy], fu2n[k, kk, iy], fv2n[k, kk, iy], fuun[k, kk, iy], fuvn[k, kk, iy], fvvn[k, kk, iy]))
        f.close()
end = time.perf_counter()
print(end-start)
