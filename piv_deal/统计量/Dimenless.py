import numpy as np
import os
# -*- coding:utf-8 -*-

m = 399
n = 65
s = m*n
path_data = "L:\\M_C\\4p\\"
para_metre = path_data + 'parameter.txt'
path_data1 = path_data + "data\\"
path_result = path_data + "Result\\Statistic\\"
path_mean = path_result + 'Velocity_mean.dat'
path_list = os.listdir(path_data1)
path_tau = path_result + 'Mu_tau.dat'
path_Reynold_Stress = path_result + 'Reynold_Stress.dat'
path_Statistic = path_result + 'Statistic.dat'
# 读取参数信息
f = open(para_metre, 'r')
line1 = f.readline()
sdy = line1.split("=")
mul = float(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
dyn = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
ps = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
pe = int(sdy[1])
f.close()
# 读取平均速度场
mean_velocity = np.loadtxt(path_mean, dtype='float', skiprows=2)
meany = mean_velocity[:, 1].reshape(m, n, order='F')[1, :]
u = mean_velocity[:, 2].reshape(m, n, order='F')
um = np.mean(u, 0)
# 读取壁面摩擦速度
u_tau = np.loadtxt(path_tau, dtype='float', skiprows=2)
U_tau = u_tau[-1]
# 读取雷诺应力场
Rs = np.loadtxt(path_Reynold_Stress, dtype='float', skiprows=2)
fuu = Rs[:, 2].reshape(m, n, order='F')
fvv = Rs[:, 3].reshape(m, n, order='F')
fuv = Rs[:, 4].reshape(m, n, order='F')
fuul = np.mean(np.sqrt(fuu)[5:m-5][:], 0) / U_tau
fvvl = np.mean(np.sqrt(fvv)[5:m-5][:], 0) / U_tau
fuvl = np.mean(fuv[5:m-5][:], 0) / U_tau / U_tau
ul = um / U_tau
yl = meany * U_tau / mul / 1000
f = open(path_Statistic, 'w')
f.write('VARIABLES= ,"<times>y<sup>+</sup></times>", "<times>u<sup>+</sup></times>", "<times>urms<sup>+</sup></times>", "<times>vrms<sup>+</sup></times>", "<times>uv<sup>+</sup></times>"\n')
f.write('Zone T=%s\t%s%d\n' % ('"1p"', ' I=', n))
for i in range(0, n):
    f.write("{}\t{}\t{}\t{}\t{}\n".format(yl[i], ul[i], fuul[i], fvvl[i], fuvl[i]))
f.close()
