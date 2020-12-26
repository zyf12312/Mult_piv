# 此程序使用对数率区公式拟合平均速度剖面，求取壁面摩擦速度
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

m = 399
n = 65
a = 0.41
b = 5.0
path_data = "L:\\M_C\\4p\\"

para_metre = path_data + 'parameter.txt'
f = open(para_metre, 'r')
line1 = f.readline()
sdy = line1.split("=")
mul = float(sdy[1])
skip = f.readline()
line1 = f.readline()
sdy = line1.split("=")
ps = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
pe = int(sdy[1])
f.close()

path_data1 = path_data + "data\\"
path_result = path_data + "Result\\Statistic\\"
path_mean = path_result + 'Velocity_mean.dat'
path_tau = path_result + 'Mu_tau.dat'
path_sta = path_result + 'information.dat'
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meany = mean_velocity[:, 1].reshape(m, n, order='F')
u = mean_velocity[:, 2].reshape(m, n, order='F')
y = meany[1, :] / 1000
um = np.mean(u[5:m - 5][:], 0)
uq = um.max()


def func(x, ue_spl):
    return mul / ue_spl * (x / ue_spl + np.exp(-a * b) * (
            np.exp(a * x / ue_spl) - 1 - a * x / ue_spl - (a * x / ue_spl) ** 2 / 2 - ((a * x / ue_spl) ** 3) / 6 - (
            a * x / ue_spl) ** 4 / 24))


popt, pcov = curve_fit(func, um[ps:pe + 1], y[ps:pe + 1], p0=0.01)
print(popt[0])
i = 0
while um[i] < 0.99 * uq:
    i = i + 1
bc = y[i - 1] + (0.99 * uq - um[i - 1]) / (um[i] - um[i - 1]) * (y[i] - y[i - 1])
ny = pe - ps + 1
e = np.zeros(ny)

for i in range(0, ny):
    e[i] = 0.00001
    ae = 1.0
    while abs(e[i] - ae) >= 1.0e-10:
        ae = e[i]
        f1 = 1 / a / e[i] + um[i + ps] / (e[i] * e[i])
        f2 = 1 / a * np.log(y[i + ps] * e[i] / mul) + b - um[i + ps] / e[i]
        e[i] = ae - f2 / f1

ue = np.mean(e)
retua = ue * bc / mul
print(ue)

f = open(path_tau, 'w')
f.write('VARIABLES= "Mu_tau(m/s)"\n')
f.write('Zone T="%s"\n' % '1p')
for i in range(0, ny):
    f.write("{}\n".format(e[i]))
f.write("{}\n".format(ue))
f.write("{}".format(popt[0]))
f.close()

f = open(path_sta, 'w')
f.write('bc=%f\n' % bc)
f.write('uq=%f\n' % uq)
f.write('ue=%f\n' % ue)
f.write('uespl=%f\n' % popt[0])
f.write('Retau=%f' % retua)
f.close()
