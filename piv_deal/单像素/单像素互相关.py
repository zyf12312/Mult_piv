# 单像素互相关，近壁，原点783像素点
# -*- coding:utf-8 -*-
import numpy as np
import time
from PIL import Image
from scipy import interpolate
import os

np.set_printoptions(suppress=True)
start = time.perf_counter()
#####################################################################################################
# 基本信息
disk = "L:\\单像素\\"
case = '1p'
para_metre = disk + case + '_parameter.txt'
f = open(para_metre, 'r')
line1 = f.readline()
sdy = line1.split("=")
t = 1 / int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
cal = float(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
xp = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
yp = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
zero = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
yp1 = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
yp2 = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
y1 = int(sdy[1])
line1 = f.readline()
sdy = line1.split("=")
y2 = int(sdy[1])
f.close()
td = 0
n = 8215
dx = 30  # 查询窗口在流向大小
dy = 6  # 查询窗口在流向大小
det_y = 1
det_x = 2
#####################################################################################################
# 文件位置
path = '%s%s\\' % (disk, case)
path_list = os.listdir(path)
path1 = disk + case + '_单像素.dat'
#####################################################################################################
# 读数据
sumI = np.zeros((yp1 - yp2, xp))
Imall = np.zeros((n, yp1 - yp2, xp))
for i in range(n):
    str1 = path + path_list[i]
    im = np.array(Image.open(str1).convert('L'), 'f')
    im = np.array(im)[yp2 - 1: yp1 - 1, :]
    Imall[i, :, :] = im
    sumI += im
print('组合完成')
#####################################################################################################
# 平均和均方根
I1a = (sumI - Imall[0, :, :]) / (n - 1)
I2a = (sumI - Imall[n - 1, :, :]) / (n - 1)
dI1 = np.zeros((yp1 - yp2, xp))
dI2 = np.zeros((yp1 - yp2, xp))
for i in range(n - 1):
    dI1 += (Imall[i, :, :] - I1a) ** 2
    dI2 += (Imall[i + 1, :, :] - I2a) ** 2
dI1 /= (n - 2)
dI2 /= (n - 2)
sigmaI1 = np.sqrt(dI1)
sigmaI2 = np.sqrt(dI2)
#####################################################################################################
# 相关计算
aa = np.zeros(y2 - y1)
c = np.zeros(int((y2 - y1) / det_y))
y = np.arange(y1, y2)
y = ((zero - y) * cal)
for i in range(y1 - yp2, y2 - yp2, det_y):
    for j in range(150, 1050, det_x):
        cor = np.zeros((dy * 2 + 1, dx + 1))
        for k in range(n - 1):
            im1 = Imall[k, i, j]
            im2 = Imall[k + 1, i - dy: i + dy + 1, j: j + dx + 1]
            cor += (im1 - I1a[i, j]) * (im2 - I2a[i - dy: i + dy + 1, j: j + dx + 1])
        R = (cor / (sigmaI1[i, j] * sigmaI2[i - dy: i + dy + 1, j: j + dx + 1])) / (n - 2)
        loc1 = np.argwhere(R == R.max())
        cy = R.max()
        yy = loc1[0][0]
        xx = loc1[0][1]
        if 1 < yy < 2 * dy - 1 and 1 < xx < dx - 1:
            c[i - (y1 - yp2)] += 1
            dxx = (np.log(R[yy, xx - 1]) - np.log(R[yy, xx + 1])) / \
                  (2 * (np.log(R[yy, xx + 1]) - 2 * np.log(R[yy, xx]) + np.log(R[yy, xx - 1])))
            rxx = xx + dxx
            aa[i - (y1 - yp2)] += rxx
            # # 二维插值
            # xi = np.linspace(-2, 2, 5)
            # yi = np.linspace(-2, 2, 5)
            # xi, yi = np.meshgrid(xi, yi)
            # zi = R[yy - 2: yy + 3, xx - 2: xx + 3]
            # inter_func = interpolate.interp2d(xi, yi, zi, kind='cubic')
            # xn = np.linspace(-1, 1, 101)
            # yn = np.linspace(-1, 1, 101)
            # Rn = inter_func(xn, yn)
            # loc2 = np.argwhere(Rn.max() == Rn)
            # rxx = xx + xn[loc2[0][1]]
            # aa[i - (y1 - yp2)] += rxx
dd = aa / c
v = (dd * cal / t) / 1000
#####################################################################################################
result = np.row_stack((y, v)).T
f = open(path1, 'w')
f.write('VARIABLES= "y(mm)", "u(m/s)" \n')
f.write('Zone T="1p" i=%d \n' % len(y))
for a in range(len(result)):
    for b in range(len(result[0])):
        f.write(str(result[a][b]))
        f.write('\t')
    f.write('\n')
f.close()
#####################################################################################################
end = time.perf_counter()
print('Running time: %s Seconds' % (end - start))
