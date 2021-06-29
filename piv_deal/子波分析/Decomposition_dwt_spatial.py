# decomposition flow stream-wise spatially by using discrete wavelet transform
import numpy as np
import time
import pywt
import csv
from Read_Parameter import read_p


###################################################################
def wrcoef(nn, coef_type, coeffs, wavename, levels):
    a, ds = coeffs[0], list(reversed(coeffs[1:]))
    if coef_type == 'a':
        return pywt.upcoef('a', a, wavename, level=levels, take=nn)
    elif coef_type == 'd':
        return pywt.upcoef('d', ds[levels - 1], wavename, level=levels, take=nn)
    else:
        raise ValueError("Invalid coefficient type:{}".format(coef_type))


###################################################################
# read parameters from file
path_data = "L:\\M_C\\1p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8  # number of parameters
nu = 1  # number of cases
bu = read_p(para_metre, p_n)
mul = float(bu[0])  # coefficient of kinematic viscosity
NN = int(bu[4])  # number of samples
m = int(bu[5])  # number of stream-wise points
n = int(bu[6])  # number of wall-normal points
##################################################################
# basic parameter
N = NN * nu  # sum of all frames
fs = 200  # frequency of sampling
dt = 1 / fs  # interval times between adjacent frames
start = time.perf_counter()
w = pywt.Wavelet('db3')
ll = pywt.dwt_max_level(m, w)+5
ss = (n, ll)
uf = np.zeros(ss)
vf = np.zeros(ss)
uft = np.zeros(ss)
vft = np.zeros(ss)
sss = np.zeros((1, m))
##################################################################
# Set path
format_fuv = '{}d'.format(m * n)
s = (m, n)
fu = np.zeros(s)
fv = np.zeros(s)
path_result = path_data + "Result\\"
path_fuv = path_result + 'Fuv\\'
path_dwt_s = path_result + 'dwt_spatial\\'
path_mean = path_result + 'Statistic\\Velocity_mean.dat'
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meanx = mean_velocity[:, 0]
meany = mean_velocity[:, 1]
for i in range(N):
    path_fu = path_fuv + "fu_{:05}.bin".format(i + 1)
    path_fv = path_fuv + "fv_{:05}.bin".format(i + 1)
    fu = np.fromfile(path_fu, dtype=format_fuv).reshape(m, n, order='F')
    fv = np.fromfile(path_fv, dtype=format_fuv).reshape(m, n, order='F')
    tu = 0.5 * np.mean(fu ** 2, axis=0)
    tv = 0.5 * np.mean(fv ** 2, axis=0)
    for j in range(n):
        uce = pywt.wavedec(fu[:, j], w, level=ll)
        vce = pywt.wavedec(fv[:, j], w, level=ll)
        for k in range(ll):
            uu = wrcoef(m, 'd', uce, w, k + 1)
            vv = wrcoef(m, 'd', vce, w, k + 1)
            uf[j, k] = 0.5 * np.mean(uu ** 2)/tu[j]
            vf[j, k] = 0.5 * np.mean(vv ** 2)/tv[j]
    uft += uf / N
    vft += vf / N
ll1 = uft.argmax(axis=1)
ll2 = vft.argmax(axis=1)
with open(path_dwt_s + 'energy_fu_5.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow(range(1, ll + 1))
    writer.writerows(uft)
with open(path_dwt_s + 'energy_fv_5.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerow(range(1, ll + 1))
    writer.writerows(vft)
path_dwt_energy_max = path_dwt_s + "energy_5.dat"
fid = open(path_dwt_energy_max, 'w')
fid.write('VARIABLES="y","uMaxscale","vMaxscale"\n')
fid.write('Zone T="1" i=%d\n' % n)
for i in range(n):
    fid.write("{}\t{}\t{}\n".format(i + 1, ll1[i], ll2[i]))
fid.close()
