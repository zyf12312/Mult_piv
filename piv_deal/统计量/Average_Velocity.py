#  此程序求取二维PIV平均速度场
import numpy as np
import os
from Read_Parameter import read_p
if __name__ == "__main__":
    read_p
path_data = "L:\\M_C\\1p\\"
para_metre = path_data + 'parameter.txt'
p_n = 8
bu = read_p(para_metre, p_n)
mul = float(bu[0])         # coefficient of kinematic viscosity
dyn = int(bu[1])           # Deviance from zero to real wall coordinate
ps = int(bu[2])            # first point in log law region
pe = int(bu[3])            # last point in log law region
NN = int(bu[4])            # number of samples
m = int(bu[5])             # number of stream-wise points
n = int(bu[6])             # number of wall-normal points
scale = float(bu[7])       # amplification coefficient

dy = -0.075*dyn

nu = 6
N = NN * nu
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

