#  此程序求取二维PIV平均速度场
import numpy as np
import os

nu = 6
NN = 8215

scale = 3.75
N = NN * nu
m = 399
n = 65
path_data = "L:\\M_C\\1p\\"
para_metre = path_data + 'parameter.txt'
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
dy = -0.075*dyn

path_data1 = path_data + "data\\"
path_result = path_data + "Result\\"
path_mean = path_result + 'Velocity_mean.dat'
path_list = os.listdir(path_data1)

path_velocity = path_data1 + path_list[0]
print(path_velocity)
V = np.loadtxt(path_velocity, dtype='float', skiprows=3, comments='DATASETAUXDATA')
m_v = np.zeros(np.shape(V))
m_v[:, 0] = V[:, 0] * scale
m_v[:, 1] = V[:, 1] * scale + dy

for i in range(0, N):
    path_velocity = path_data1 + path_list[i]
    V = np.loadtxt(path_velocity, dtype='float', skiprows=3, comments='DATASETAUXDATA')
    m_v[:, 2:4] = m_v[:, 2:4] + V[:, 2:4]
m_v[:, 2:4] = m_v[:, 2:4] / N * scale

f = open(path_mean, 'w')
f.write('VARIABLES= "x(mm)", "y(mm)", "u(m/s)", "v(m/s)" \n')
f.write('Zone T="%d" i=%d j=%d\n' % (nu, m, n))
for a in range(len(m_v)):
    for b in range(len(m_v[0])):

        f.write("{}".format(m_v[a][b]))
        f.write('\t')

    f.write('\n')

f.close()
